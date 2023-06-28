"""
Algorithm for filtering low SNR scintillating sources

"""


import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
#import time as t
from windows import tukey
#import scipy.stats as st
from utils import maxi_sub,maxi_unsub,trackcoord,ang_dist_mult
from scipy.signal import find_peaks, peak_widths, peak_prominences

vlss=np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt")
def search_vlss(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist_mult(ra_in, dec_in, vlss[:,0], vlss[:,1])
    par = np.where((vlss[:,2] > 10.0) & (diff < 3.0))[0]

    return len(par)



def filter_snr(inf):

    num1 = 0
    stokes = [0]
    list1 = np.zeros(inf.shape)
    for hour in range(0,24):
        #find the corresponding .hdf5 file for the time 
        path1 = '/leo/savin/ovro-lwa/hdf5_files/aug10/'+'ovro_10aug_01band_'+str(format(hour,'02d'))+'h.hdf5'
        files1  = sorted(glob.glob(path1))
        hour_index = np.where((inf[:,3] == hour))[0]
        if len(hour_index) == 0:
           print "no transients at utc" +str(hour)
           continue

        for filename1 in files1:
            print filename1
            hf = h5py.File(filename1,'r')
            img = hf.get('image')
            tm = hf.get('time')
            ints = img.shape[0] #number of integrations
            ngrid  = 4096 #size of images (x-axis)
            psize = 1.875/60.0 #angular size of a pixel (at zenith)
            t_comp = tm[:,1]*3600.0 + tm[:,2]*60.0 + tm[:,3]
            

            for i in hour_index:
                ra = inf[i,0]
                dec = inf[i,1]
                MJD = inf[i,2]
                h = inf[i,3]
                m = inf[i,4]
                s = inf[i,5]
       
                mint = int(m)
                sec = s+(m-mint)*60.0
                if sec >= 60.0:
                   sec = int(sec%60.0)
                   mint += 1
                   if mint >= 60:
                      mint = mint%60
                      h += 1 
                sec = int(sec)
                h  = int(h) 
                t_in = h*3600.0 + mint*60.0 + sec      
                peak = np.where((abs(t_comp-t_in) < 3))[0]
                if len(peak)==0:
                   #print"requested time not in dataset"
                   print inf[i,:]
                   continue 
                if len(peak) > 1:
                   peak = int(np.median(peak))
                    
                #signal = np.arange(max(0,peak-305),min(peak+305,ints))
                signal = np.arange(ints)
                noise = np.where((signal < peak-5)| (signal > peak+5))[0]
                event = np.where((signal > peak-3)&(signal < peak+3))[0]
    
                intg = len(signal)
                if intg == 0:
                   continue
                light=np.zeros((intg),dtype=float)
                data_filt = np.zeros((intg),dtype=np.complex64)
       
                lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
                lon_ovro = -118.28166  #Station Latitude (of LWA-ovro) 

                for ln,j in enumerate(signal):
                    time = tm[j,:]
                    xpix,ypix,altpix=trackcoord(np.array([ra]),np.array([dec]),time[0],time[1],time[2],time[3],ngrid,psize,lat_ovro,lon_ovro)
                    data_i = img[j,:,:]
                    xa, ya = maxi_unsub(data_i,xpix,ypix)
                    light[ln] = data_i[xa,ya] 
            
               
                y1 = np.fft.fft(light)
                y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
                data_filt = np.fft.ifft(y1)
             
                 
                sig_noise = np.std(data_filt.real[noise])
                mean = np.mean(data_filt.real[noise])
                   
                snr = (data_filt.real[event].max()-mean)/sig_noise               
                
                """ 
                if (snr > 5.0):
                   list1[num1,:] = inf[i,:]
                   num1 += 1
                
                """
                dat_noise = data_filt.real[noise] - mean
                sig  = np.std(dat_noise)

                peaks,_ = find_peaks(dat_noise, height = 3.5*sig, distance = 3)
                top = len(peaks)
                if (snr > 5.0):
                   if top < 2:
                      list1[num1,:] = inf[i,:]
                      num1 += 1
                   else:
                       if (search_vlss(inf[i,0],inf[i,1]) == 0) :
                          list1[num1,:] = inf[i,:]
                          num1 += 1

                  

    
            hf.close() 

    return list1[:num1,:]    

if __name__ == '__main__':
   
   mjd = 58000
   path_file = 'transients/'+str(mjd)+'/*.txt'
   files_txt = sorted(glob.glob(path_file))
   for filetxt in files_txt:
       print filetxt
       inp = np.loadtxt(filetxt)
       out  =  filter_snr(inp)
       np.savetxt('snr_transients/'+str(mjd)+'/'+filetxt, out, '%10.3f')
