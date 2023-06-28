"""
Filtering low snr scintillation events
"""

import numpy as np
import sys,os
#import matplotlib.pyplot as plt
import glob
import h5py
#import time as t
from windows import tukey
#import scipy.stats as st
from utils_new import maxi_sub,maxi_unsub,trackcoord,ang_dist_mult
from scipy.signal import find_peaks, peak_widths, peak_prominences

vlss=np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt")
def search_vlss(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist_mult(ra_in, dec_in, vlss[:,0], vlss[:,1])
    par = np.where((vlss[:,2] > 10.0) & (diff < 3.0))[0]

    return len(par)



def filter_snr(inf,mjd):

    num1 = 0
    stokes = [0]
    list1 = np.zeros(inf.shape)
    for hour in range(0,24):
        #find the corresponding .hdf5 file for the time 
        path1 = '/leo/savin/ovro-lwa/hdf5_files/aug16/'+'ovro_16aug_01band_'+str(format(hour,'02d'))+'h.hdf5'
        files1  = sorted(glob.glob(path1))
        hour_index = np.where((inf[:,2] == mjd) & (inf[:,3] == hour))[0]
        print len(hour_index)
        if len(hour_index) == 0:
           print "no transients at utc" +str(hour)
           continue
        peak_ar = np.zeros(len(hour_index))
        for filename1 in files1:
            print filename1
            try:
               hf = h5py.File(filename1,'r')
               img = hf.get('image')
               tm = hf.get('time')
            except IOError:
               continue      
            ints = img.shape[0] #number of integrations
            ngrid  = 4096 #size of images (x-axis)
            psize = 1.875/60.0 #angular size of a pixel (at zenith)
            lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
            lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)
            t_comp = tm[:,1] + tm[:,2]/60.0 + tm[:,3]/3600.00
            

            for bi,i in enumerate(hour_index):
                #print 'finding peaks'
                #print inf[i,:]
                #ra = inf[i,0]
                #dec = inf[i,1]
                #MJD = inf[i,2]
                h = inf[i,3]
                m = inf[i,4]
                s = inf[i,5]
                t_in = h + m/60.0 + s/3600.00      
                pk = np.where((abs(t_comp-t_in)*3600.00 < 20))[0]
                #print pk
                if len(pk) > 0:
                   peak_ar[bi] = int(np.median(pk))
                else:
                   peak_ar[bi] = 1000 

            coord_ar = np.zeros((len(hour_index),ints,2),dtype = int)
            for m,j in enumerate(hour_index):
                #print 'finding coordinates'
                ra = inf[j,0]
                dec = inf[j,1]
                #print ra,dec
                ra_ar = np.ones(ints)*ra
                dec_ar = np.ones(ints)*dec
                xpix,ypix,altpix=trackcoord(ra_ar, dec_ar, tm[:,0], tm[:,1], tm[:,2], tm[:,3], ngrid, psize, lat_ovro, lon_ovro)
                #print xpix,ypix
                #print np.concatenate((xpix,ypix),axis=1).shape
                coord_ar[m,:,0] = xpix#np.concatenate((xpix,ypix),axis=1)
                coord_ar[m,:,1] = ypix

            light_ar = np.zeros((len(hour_index),ints)) 
            for k in range(ints):
                #print 'collecting data'
                #print k
                data = img[k,:,:]
                for n in range(len(hour_index)):
                    xa,ya = coord_ar[n,k,:]
                    if (0 <= xa < 4096) & (0 <= ya < 4096):
                       #print xa,ya
                       light_ar[n,k] = data[xa,ya] 

            for l,q in enumerate(hour_index):
                #print 'filtering'
                #print l
                peak = peak_ar[l]
                if peak == 1000:
                   print 'No peak'
                   continue
                signal = np.arange(ints)
                noise = np.where((signal < peak-5)| (signal > peak+5))[0]
                event = np.where((signal > peak-3)&(signal < peak+3))[0]
                
                light = light_ar[l,:]
                
                y1 = np.fft.fft(light)
                y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
                data_filt = np.fft.ifft(y1)
                #plt.plot(signal,data_filt)
                #plt.show()
                 
                sig_noise = np.std(data_filt.real[noise])
                mean = np.mean(data_filt.real[noise])
                   
                snr = (data_filt.real[event].max()-mean)/sig_noise               
                
                #print snr
                dat_noise = data_filt.real[noise] - mean
                sig  = np.std(dat_noise)

                peaks,_ = find_peaks(dat_noise, height = 3.5*sig, distance = 3)
                top = len(peaks)
                if (snr > 5.0):
                   if top < 2:
                      list1[num1,:] = inf[q,:]
                      num1 += 1
                   else:
                       if (search_vlss(inf[q,0],inf[q,1]) == 0) :
                          list1[num1,:] = inf[q,:]
                          num1 += 1

                  

    
            hf.close() 

    return list1[:num1,:]    

if __name__ == '__main__':
   
   mjd = 58346
   path_file = 'transients/'+str(mjd)+'/*.txt'
   files_txt = sorted(glob.glob(path_file))
   for filetxt in files_txt:
       print filetxt
       basetxt = os.path.basename(filetxt)
       #print basetxt
       inp = np.loadtxt(filetxt)
       out  =  filter_snr(inp,mjd)
       if out.shape[0] != 0:
          np.savetxt('snr_transients/'+str(mjd)+'/'+basetxt, out, '%10.3f')
