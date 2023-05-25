import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import h5py
from utils import trackcoord,maxi_unsub
from scipy.stats import scoreatpercentile as sc
import matplotlib.animation as animation
from sliding_rfi_flagger import main
from utils import mask_image,maxi_sub,maxi_unsub,trackcoord

ra, dec, mjd, h, m, s, filename = sys.argv[1:]

ra = float(ra)
dec = float(dec)
mjd = int(float(mjd))
h = int(float(h))
m = int(float(m))
s = int(float(s))
#stokes = int(stokes)

hf=h5py.File(filename,'r')
dset1 = hf.get('image')
dset2 = hf.get('time')
header = hf.get('header')
chan = hf.get('channel')[0,:]
ints = dset1.shape[0] #number of integrations
ngrid  = header[0,1] #size of images (x-axis)
psize = header[0,0] #angular size of a pixel (at zenith)

print(chan.shape)

lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)





time_in = h*3600.00 + m*60.0 + s

time_dat = dset2[:,1]*3600.00 + dset2[:,2]*60.0 + dset2[:,3]
peak = np.where((np.abs(time_dat - time_in) < 13.0))[0]
print (peak)

peak_ind = int(np.median(peak))

noise = np.array(range(max(0,peak_ind-10),max(0,peak_ind-6))+range(min(ints,peak_ind+6),min(ints,peak_ind+10)))

noise_im  = np.mean(dset1[noise,:,:,:],axis=0)


if len(peak) == 0:
   sys.exit(" File grabbed is not right")
      

def on_press(event):
    if event.key.isspace():
        if anim.running:
            anim.event_source.stop()
        else:
            anim.event_source.start()
        anim.running ^= True


#making movie 
fig = plt.figure()

listims = []
fig.canvas.mpl_connect('key_press_event', on_press)

for index in peak :
    #print index
    time = dset2[index,:]

    xp,yp,alt = trackcoord(np.array([ra]),np.array([dec]),time[0],time[1],time[2],time[3],ngrid,psize,lat_ovro,lon_ovro)
    
    print(" xp, yp : [ %d, %d] " % (xp,yp))

    dat_img = dset1[index,:,:,:]
    for chan_ind in range(dat_img.shape[-1]): 
        im_chan = dat_img[:,:, chan_ind] - noise_im[:,:, chan_ind]
        vmin = sc(im_chan,0.1)
        vmax = sc(im_chan,99.9)
        im = plt.pcolormesh(im_chan.T,cmap ='jet',vmin = vmin,vmax=vmax)
        str1 = "loc:(%d,%d)" % (xp,yp) 
        str2 = "%0.3f MHz" % (chan[chan_ind])
        freq = plt.text(200,240,str2)
        loc = plt.text(xp,yp+40,str1)
        #circ = plt.scatter(xp, yp, 1000, facecolors = None)
        #plt.legend()
        #plt.title(" %f MHz" % (chan[chan_ind]))

        listims.append([im, freq, loc])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
anim.save('mra5_spectra_lowres.mp4')
anim.running = True
plt.show()

