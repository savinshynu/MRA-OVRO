import numpy as np
import sys,glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import os
from scipy.stats import scoreatpercentile as sp

path = 'events/mra1/nearcor/*.hdf5'

files = sorted(glob.glob(path))
capt = 13
outfile = h5py.File('events/mra1/nearcor/ovro_aug10_mra1.hdf5','a')
# Initializing a hdf5  dataset
dset = outfile.create_dataset('image',(capt,4096,4096,4),dtype='float32') #just stokes I
tset = outfile.create_dataset('time',(capt,4),dtype='float32') #just 
hset = outfile.create_dataset('header',(1,2),dtype='float32') #just
#chset = outfile.create_dataset('channel',(1,436),dtype='float32') #just


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
    #start_freq, stop_freq, bandwidth = head[:3]
    #nchan = 109
    #freq = np.linspace(start_freq,stop_freq,nchan)
    #print tset1[:,1]
    if i == 0:
       tset[:,:] = tset1
       hset[0,:] = np.array([psize,ngrid])
    dset[:,:,:,i] = dset1
    print tset1[:,3]

outfile.close()
