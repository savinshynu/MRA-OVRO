import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
from windows import tukey
from utils import maxi_unsub, trackcoord, avgp_hres, sub_res_peak
from scipy.signal import find_peaks, peak_widths, peak_prominences
from sliding_rfi_flagger import main


filename = sys.argv[1]
print filename
print 'getting data from hdf5 files'

hf=h5py.File(filename,'r')
dset1 = hf.get('image')
dset2 = hf.get('time')
dset3 = hf.get('header')
hdr  = dset3[0,:]

ints =  dset1.shape[0]

lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)


ra = np.array([float(sys.argv[2])])
dec = np.array([float(sys.argv[3])])
mjd=float(sys.argv[4])
h=float(sys.argv[5])
m=float(sys.argv[6])
s=float(sys.argv[7])
#stokes= int(sys.argv[4])

#size = hdr[4] 
#psize = hdr[3]
size  = dset3[0,1] #size of images (x-axis)
psize = dset3[0,0] #angular size of a pixel (at zenith)

print ints,size, psize

t_comp = dset2[:,1] + (dset2[:,2]/60.0) + (dset2[:,3]/3600.00)
t_in = h + (m/60.0) + (s/3600.00)        
        
try:
    peak = np.where((abs(t_comp-t_in)*3600.00 < 6))[0]
except IndexError:
    sys.exit('No peak value')
if len(peak) == 0: 
   sys.exit('No peak time') 

peak_ind = int(np.median(peak))
signal = np.arange(max(0,peak_ind-50),min(ints,peak_ind+50))



k=0
light=np.zeros(len(signal))

for n,j in enumerate(signal):

    data = dset1[j,:,:,:]
    time = dset2[j,:]
    xpix,ypix,altpix=trackcoord(ra, dec, time[0], time[1], time[2], time[3], size, psize, lat_ovro, lon_ovro)
    
    dat_img = np.mean(data[:,:,:],axis=2)
    xa, ya = maxi_unsub(dat_img,xpix,ypix)
    
    gd = main(data[xa,ya,:])
    #gd = list(set(gd)-set(range(218,248)))
    dat_img_flag = np.mean(data[:,:,gd],axis=2)


    #light[j] =  dat_img_flag[xa,ya] #avgp(xa,ya,data)
    light[n] =  sub_res_peak(xa,ya,dat_img_flag)
    print n, light[n],time[1],time[2],time[3],xpix,ypix,xa,ya,(altpix*180.0/np.pi)
    


#y1 = np.fft.fft(light[:k,1])
y1 = np.fft.fft(light)
y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
data_filt = np.fft.ifft(y1)

#save the light curves
np.savetxt("lc_mra5_subt.txt", light)

# plotttingpart 

fig,ax= plt.subplots()
ax.plot(range(light.shape[0]),light,color='red', marker='.', linestyle='solid')
#ax.plot(range(data_filt.shape[0]),data_filt,color='blue', marker='.', linestyle='solid')
ax.set(xlabel= 'time (integrations)',ylabel= 'Arbitrary Power',title='Transient light curve')
ax.grid()
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
#print data_filt.real[242-2:242+2].max()/np.std(data_filt[range(242-55,242-5)+range(242+5,242+55)])
#print data_filt.real[242-2:242+2].max()/np.std(data_filt[:])
"""


