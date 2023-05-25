import numpy as np
import sys,glob
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#import pyfits
from matplotlib.figure import Figure
import h5py
import os
#import time
from scipy.stats import scoreatpercentile as sp


path1 = '/leo/savin/ovro-lwa/hdf5_files/aug10/'+'ovro_10aug_01band_'+'*.hdf5'
files1  = sorted(glob.glob(path1))

for filename in files1:
    print filename
    #print 'getting data from hdf5 files' 
    hf=h5py.File(filename,'r')
    dset=hf.get('image')
    
    hour = os.path.basename(filename).split('_')[3][:2]
    path2 = '/leo/savin/ovro-lwa/images_new/%s' % (hour)
    try:
        os.mkdir(path2)
    except OSError:
        print 'directry already exists'
    #making movie in a different way
    #os.system('mkdir images/10/')
    #stokes = 0
    for i in range(dset.shape[0]):
        im_dat = np.transpose(dset[i,:,:])
        vmin = sp(im_dat,0.5)
        vmax = sp(im_dat,99.5)
        fname = path2+'/fig%03d.png' %(i)
        #plt.title('%d' % i)
        plt.imsave(fname,im_dat,vmin = vmin, vmax = vmax, cmap = 'jet', origin = 'lower')
        #plt.pcolormesh(np.transpose(im_dat),cmap='jet',vmin = vmin,vmax =vmax) 
        #plt.savefig('images/fig%03d.png' %(i))
        print i

    #os.system('ffmpeg -framerate 5 -i images_new/10/fig%03d.png noiter_ovro_20180810_10_i.mp4')
