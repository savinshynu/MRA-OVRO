"""
Plot light curves by smoothening the data

"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
from windows import tukey
from utils_new import maxi_unsub, trackcoord
from scipy.signal import find_peaks, peak_widths, peak_prominences

filename = sys.argv[1]
print filename
print 'getting data from hdf5 files'

try:
   hf=h5py.File(filename,'r')
   dset1 = hf.get('image')
   time = hf.get('time')
except IOError:
   print "bad file" 
#dset3 = hf.get('header')
#hdr  = dset3[0,:]
intg =  dset1.shape[0]
xSize = 4096
ySize = 4096   

lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)


ra = np.array([float(sys.argv[2])])
dec = np.array([float(sys.argv[3])])
#MJD=float(sys.argv[4])
#h=float(sys.argv[5])
#m=float(sys.argv[6])
#s=float(sys.argv[7])
#stokes= int(sys.argv[4])

size = 4096 
psize = 1.875/60.0 

ra_ar = np.ones(intg)*ra
dec_ar = np.ones(intg)*dec

xpix,ypix,altpix=trackcoord(ra_ar, dec_ar, time[:,0], time[:,1], time[:,2], time[:,3], size, psize, lat_ovro, lon_ovro)
#print xpix
#print ypix
it = np.arange(intg)
light  = np.zeros(intg)
for i in it:
    light[i] = dset1[i,xpix[i],ypix[i]]
    print i, xpix[i],ypix[i]
    
y1 = np.fft.fft(light)
y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
data_filt = np.fft.ifft(y1)

plt.plot(it,data_filt)
plt.show()

"""
dat_noise = data_filt.real - np.mean(data_filt.real)
sig  = np.std(dat_noise)
#print sig
#print np.mean(dat_noise)

peaks,_ = find_peaks(dat_noise, height = 3.0*sig, distance = 3)
print peaks.shape
#prom = peak_prominences(dat_noise, peaks, wlen = 37 )
widths = np.zeros((len(peaks),4))
for i,pk in enumerate(peaks):
    width = peak_widths(dat_noise, np.array([pk]), rel_height= 1-((3.0*sig)/dat_noise[pk]))
    print width
    widths[i,:] = width
#width = peak_widths(dat_noise, peaks, rel_height = 0.85 )

#print widths

plt.hlines(widths[:,1],widths[:,2],widths[:,3], color="C2")
plt.plot(range(len(dat_noise)),dat_noise)
plt.scatter(peaks,dat_noise[peaks])
plt.show()
"""

