import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import h5py
from utils import trackcoord,maxi_unsub
from scipy.stats import scoreatpercentile as sc
import matplotlib.animation as animation

filename = sys.argv[1]

hf=h5py.File(filename,'r')
dset1 = hf.get('image')
dset2 = hf.get('time')
#dset3 = hf.get('header')

#hd = dset3[0,:]
ints =  dset1.shape[0]

#im  = np.transpose(np.mean(dset1[1,:,:,:],axis=2))
#plt.pcolormesh(im,cmap ='jet')
#plt.show()

def on_press(event):
    if event.key.isspace():
        if anim.running:
            anim.event_source.stop()
        else:
            anim.event_source.start()
        anim.running ^= True


#making movie 
fig = plt.figure(figsize=(8,8))

listims = []
fig.canvas.mpl_connect('key_press_event', on_press)

for index in range(ints) :
    #print index
    im_int  = np.transpose(np.mean(dset1[index,:,:,:],axis=2))
    
    vmin = sc(im_int,0)
    vmax = sc(im_int,100)
    im = plt.pcolormesh(im_int,cmap ='jet',vmin = vmin,vmax=vmax)
    #im = ax2.imshow(im_int,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
anim.save('/home/savin/plots_ovro/mra1/mra1_highres.mp4')
anim.running = True
plt.show()

