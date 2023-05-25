import glob
import os
import numpy as np
import h5py
import time
import sys


day = 'aug15'

path3 = '/leo/savin/ovro-lwa/low_res/avg_hdf5/'+day+'/'
try: 
    os.mkdir(path3)
except OSError:
               print 'directry already exists'

   
path1 = '/leo/savin/ovro-lwa/low_res/hdf5_files/'+day+'/*.hdf5'
files1 = sorted(glob.glob(path1))

for filename in files1:
   
    print filename
    #print 'getting data from hdf5 files'
    hf = h5py.File(filename,'r')
    imset = hf.get('image')
    tset = hf.get('time')
    hset = hf.get('header')
    head = hset[0,:]

    ints = imset.shape[0] #number of integration

    ngrid = head[4] #size of images (x-axis)
    psize = head[3] #angular size of a pixel (at zenith)
   
    file_ext = os.path.basename(filename)
    h5_ext = os.path.splitext(file_ext)[0]
  
    outfile = h5py.File(path3+h5_ext+'.hdf5','a') 
    dset1 = outfile.create_dataset('image',(ints,ngrid,ngrid),dtype='float32')
    dset2 = outfile.create_dataset('time',(ints,4),dtype='float32')
    dset3 = outfile.create_dataset('header',(ints,5),dtype='float32')
   
    for i in range(ints):
    
        # Copying the array from main hdf5 to avg dataset
        dset1[i,:,:] = np.mean(imset[i,:,:,:],axis=2)
        dset2[i,:] = tset[i,:]
        dset3[i,:] = hset[i,:]
      
        
         
    # closing the write out file
    outfile.close()
    
