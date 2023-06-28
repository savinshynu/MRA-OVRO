"""
View the peak image of an event

"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
from scipy.stats import scoreatpercentile as sc
from utils import maxi_sub,maxi_unsub,trackcoord,ang_dist_mult,mask_image

ra = float(sys.argv[1])
dec = float(sys.argv[2])
mjd = float(sys.argv[3])
h = float(sys.argv[4])
m = float(sys.argv[5])
s = float(sys.argv[6])

capt = float(sys.argv[7])

if capt not in [13,26,65]:
   sys.exit("Wrong integration length")

if capt == 13:
   diff_tm = 6
elif capt == 26:
   diff_tm = 14
elif capt == 65:
   diff_tm = 33


filename = sys.argv[8]
print filename
hf = h5py.File(filename,'r')
img = hf.get('image')
tm = hf.get('time')
ints = img.shape[0] #number of integrations
ngrid  = 4096 #size of images (x-axis)
psize = 1.875/60.0 #angular size of a pixel (at zenith)

lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
lon_ovro = -118.28166  #Station Latitude (of LWA-ovro

t_comp = tm[:,1] + tm[:,2]/60.0 + tm[:,3]/3600.00
            
t_in = h + m/60.0 + s/3600.00      

#diff_tm =60
peak = np.where((abs(t_comp-t_in)*3600.00 < diff_tm))[0]
print peak                
if len(peak)==0:
   sys.exit("requested time not in dataset")
                
if len(peak) > 0:
   peak_ind = int(np.median(peak))

noise = np.array(range(max(0,peak_ind-6),max(0,peak_ind-3))+range(min(ints,peak_ind+3),min(ints,peak_ind+6)))
#noise = np.array(range(max(0,peak_ind-10),max(0,peak_ind-6))+range(min(ints,peak_ind+6),min(ints,peak_ind+10)))
#sub_dat = np.mean(img[peak,:,:],axis=0) - np.mean(img[noise,:,:],axis=0)#np.mean(img[peak_ind-2:peak_ind,:,:],axis=0)
sub_dat = np.mean(img[peak,:,:],axis=0) - np.mean(img[peak_ind-2:peak_ind,:,:],axis=0)

im_mask, im_mask_nohrz = mask_image(sub_dat,tm[peak_ind,:],psize,lat_ovro,lon_ovro)

vmin = sc(im_mask_nohrz,0)
vmax = sc(im_mask_nohrz,99.99)

plt.pcolormesh(np.transpose(im_mask_nohrz),vmin=vmin,vmax =vmax,cmap='jet')
plt.colorbar()
plt.show()
                  

    
hf.close() 

