import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import pyfits
import h5py
import os
#import time
from scipy.stats import scoreatpercentile as sp
from utils import trackcoord, getcoord, gethrz, mask_image 

filename = sys.argv[1]

print filename
print 'getting data from hdf5 files' 
hf = h5py.File(filename,'r')
dset = hf.get('image')
tset = hf.get('time')
print dset.shape
im = np.mean(dset[198:199,:,:],axis = 0) #- np.mean(dset[176:178,:,:],axis = 0)
vmin = sp(im,1)
vmax = sp(im,99.5)
plt.pcolormesh(np.transpose(im),cmap='jet',vmin = vmin,vmax =vmax)
plt.show()
#print im[1785,3665] #,im.max(),im.min()
ra = np.array([299.8])
dec = np.array([40.8])
#ra, dec = [350.8, 58.8]
#time = np.array([58340, 02, 44, 40])
time = tset[178,:]

psize = 1.875/60.0
lat_ovro = 37.055166
lon_ovro =  -118.28116
size = 4096

x, y, altpix  = trackcoord(ra,dec,time[0],time[1],time[2],time[3],size,psize,lat_ovro,lon_ovro)

print x,y, altpix
t1 = np.ones(2)*time[0]
t2 = np.ones(2)*time[1]
t3 = np.ones(2)*time[2]
t4 = np.ones(2)*time[3]
ra1, dec1 = getcoord(np.array([872.45553,1000]),np.array([2559.471343,1000]),size,psize,t1,t2,t3,t4,lat_ovro,lon_ovro) 

print ra1, dec1


x1 = np.arange(0,4096)
y1 = np.ones(4096)*2048
mjd = np.ones(4096)*58340
h = np.ones(4096)*2
m = np.ones(4096)*44
s = np.ones(4096)*40

#ra1, dec1 = getcoord(x1,y1,size,psize, mjd, h, m, s, lat_ovro, lon_ovro)
#print ra1, dec1
az, alt = gethrz(x1,y1,size,psize, mjd, h, m, s, lat_ovro, lon_ovro)

#plt.plot(range(4096),alt)
#plt.show()

im_mask_hrz, im_mask_nohrz = mask_image(im,time,psize,lat_ovro,lon_ovro)

#plt.pcolormesh(np.transpose(im),cmap='jet',vmin = vmin,vmax =vmax)
#plt.show()
vmin1 = sp(im_mask_nohrz,0.1)
vmax1 = sp(im_mask_nohrz,99.9)

plt.pcolormesh(np.transpose(im_mask_nohrz),cmap='jet',vmin = vmin,vmax =vmax)
plt.show()
"""
#making movie 
stokes = 0
ims = []
fig = plt.figure()
for i in range(170,190):
    im_dat = dset[i,:,:,stokes]
    vmin = sp(im_dat,1)
    vmax = sp(im_dat,99.5)
    plt.title("%d" % i)
    im_plot = plt.pcolormesh(np.transpose(im_dat),cmap='jet',vmin = vmin,vmax =vmax)
    #plt.title("%d" % i)
    ims.append([im_plot])
    #del im_plot
    print i #sys.getsizeof(im_plot)/1e+6
    #fig1.clear()
    #plt.close(fig1)
  
ani = animation.ArtistAnimation(fig, ims, interval=500, repeat_delay=1000)
ani.save('event1.avi')
plt.show()
"""

"""

#making movies
def animate(frame):
    im_plot = None
    print frame
    stokes = 0
    
    im_dat = dset[frame,:,:,stokes]
    vmin = sp(im_dat,1)
    vmax = sp(im_dat,99.5)
    im_plot = plt.pcolormesh(np.transpose(im_dat),cmap='jet',vmin = vmin,vmax =vmax) 
    return im_plot

fig = plt.figure()
ani = animation.FuncAnimation(fig, animate,np.arange(dset.shape[0]), interval=40, repeat_delay=1000)
ani.save('ovro_10aug2018_01band_00h_i.mp4')


#making movie 
part = np.arange(0,dset.shape[0],70)
stokes = 0
for val in part:
    print val
    fig = plt.figure()
    ims = []
    for i in range(val,min(val+70,dset.shape[0])):
        im_dat = dset[i,:,:,stokes]
        vmin = sp(im_dat,1)
        vmax = sp(im_dat,99.5)
        im_plot = plt.pcolormesh(np.transpose(im_dat),cmap='jet',vmin = vmin,vmax =vmax) 
        ims.append([im_plot])
        plt.clf()
        print i

    ani = animation.ArtistAnimation(fig, ims, interval=40, blit = True, repeat_delay=1000)
    ani.save('ovro_10aug2018_01band_00h_i_part'+str(int(val/70))+'.mp4')
    ims.clear()
    #plt.show()
    fig.clf()




#making movie in a different way
os.system('mkdir images')
stokes = 0
for i in range(dset.shape[0]):
    im_dat = dset[i,:,:,stokes]
    vmin = sp(im_dat,1)
    vmax = sp(im_dat,99.5)
    plt.pcolormesh(np.transpose(im_dat),cmap='jet',vmin = vmin,vmax =vmax) 
    plt.savefig('images/fig%03d.png' %(i))
    print i
"""
