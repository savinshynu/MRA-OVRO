"""
Combine multiple low resolution images in the low band

"""

import numpy as np
import sys,glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import os
from scipy.stats import scoreatpercentile as sp

path = '/leo/savin/ovro-lwa/low_res/events/mra6/lowres/*.hdf5'

files = sorted(glob.glob(path))
#print files
hfi = h5py.File(files[0],'r')
ti = hfi.get('time')
capt = ti.shape[0]

outfile = h5py.File('/leo/savin/ovro-lwa/low_res/events/mra6/lowres/ovro_aug16_mra6.hdf5','a')
# Initializing a hdf5  dataset
dset = outfile.create_dataset('image',(capt,256,256,436),dtype='float32') #just stokes I
tset = outfile.create_dataset('time',(capt,4),dtype='float32') #just 
hset = outfile.create_dataset('header',(1,2),dtype='float32') #just
chset = outfile.create_dataset('channel',(1,436),dtype='float32') #just


for i,filename in enumerate(files):

    print filename
    hf = h5py.File(filename,'r')
    dset1 = hf.get('image')
    tset1 = hf.get('time')
    header = hf.get('header')
    head  =  header[0,:]
    ints = dset1.shape[0] #number of integrations
    ngrid  = head[4] #size of images (x-axis)
    psize = head[3] #angular size of a pixel (at zenith)
    start_freq, stop_freq, bandwidth = head[:3]
    nchan = 109
    freq = np.linspace(start_freq,stop_freq,nchan)
    print tset1[:,1]
    if i == 0:
       tset[:,:] = tset1
       hset[0,:] = np.array([psize,ngrid])
    dset[:,:,:,i*109:(i+1)*109] = dset1
    chset[0,i*109:(i+1)*109] = freq
    print tset1[:,3]

outfile.close()
