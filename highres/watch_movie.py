"""
Make a movie from high resolution images
"""

import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import h5py
from utils import trackcoord,maxi_unsub
from scipy.stats import scoreatpercentile as sc
import matplotlib.animation as animation

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

intg =  dset1.shape[0]
xSize = 4096
ySize = 4096

lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)

size = 4096
psize = 1.875/60.0




time_in = h*3600.00 + m*60.0 + s
dur  = 30.0
time_dat = dset2[:,1]*3600.00 + dset2[:,2]*60.0 + dset2[:,3]
peak = np.where((time_dat < time_in + dur) & (time_dat > time_in - dur))[0]

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
    print index
    im_int  = np.transpose(dset1[index,:,:])
    
    vmin = sc(im_int,0)
    vmax = sc(im_int,100)
    im = plt.pcolormesh(im_int,cmap ='jet',vmin = vmin,vmax=vmax)
    #im = ax2.imshow(im_int,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
#ani.save('movie_source.mp4')
anim.running = True
plt.show()

