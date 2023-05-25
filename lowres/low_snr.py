import numpy as np
import sys,os
import matplotlib.pyplot as plt
import glob
import h5py
#import time as t
from windows import tukey
#import scipy.stats as st
from utils import avgp,maxi_sub,maxi_unsub,trackcoord,ang_dist_mult
from scipy.signal import find_peaks, peak_widths, peak_prominences

vlss=np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt")
def search_vlss(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist_mult(ra_in, dec_in, vlss[:,0], vlss[:,1])
    par = np.where((vlss[:,2] > 10.0) & (diff < 3.0))[0]

    return len(par)


def search_vlss_ini(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist_mult(ra_in, dec_in, vlss[:,0], vlss[:,1])
    #print diff
    par = np.where((vlss[:,2] > 50.0) & (diff < 3.0))[0]
    #print par 
    return len(par)

def filter_snr(inf,mjd_day):

    num1 = 0
    stokes = [0]
    list1 = np.zeros((inf.shape[0],inf.shape[1]+1))
    for hour in range(0,24):
        #find the corresponding .hdf5 file for the time 
        path1 = '/leo/savin/ovro-lwa/low_res/avg_hdf5/'+mjd_day+'/ovro_'+mjd_day+'_01band_'+str(format(hour,'02d'))+'h.hdf5'
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
            hd = hf.get('header')
            head = hd[0,:]
            ints = img.shape[0] #number of integrations
            ngrid  =  head[4] #size of images (x-axis)
            psize = head[3] #angular size of a pixel (at zenith)
            t_comp = tm[:,1] + tm[:,2]/60.0 + tm[:,3]/3600.00
            

            for i in hour_index:
                ra = inf[i,0]
                dec = inf[i,1]
                MJD = inf[i,2]
                h = inf[i,3]
                m = inf[i,4]
                s = inf[i,5]
       
                t_in = h + m/60.0 + s/3600.00   
                peak = np.where((abs(t_comp-t_in)*3600 < 13))[0]
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
                #data_filt = np.zeros((intg),dtype=np.complex64)
       
                lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
                lon_ovro = -118.28166  #Station Latitude (of LWA-ovro) 

                for ln,j in enumerate(signal):
                    time = tm[j,:]
                    xpix,ypix,altpix=trackcoord(np.array([ra]),np.array([dec]),time[0],time[1],time[2],time[3],ngrid,psize,lat_ovro,lon_ovro)
                    data_i = img[j,:,:]
                    xa, ya = maxi_unsub(data_i,xpix,ypix)
                    light[ln] = avgp(xa,ya,data_i) #data_i[xa,ya] 
            
               
                y1 = np.fft.fft(light)
                y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
                data_filt = np.fft.ifft(y1)
             
                 
                sig_noise = np.std(data_filt.real[noise])
                mean = np.mean(data_filt.real[noise])
                   
                snr = (data_filt.real[event].max()-mean)/sig_noise               
                
                #print snr
                
                if (snr > 5.0):
                   list1[num1,:6] = inf[i,:]
                   list1[num1,6] = snr
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
                
                """  

    
            hf.close() 

    return list1[:num1,:]    

if __name__ == '__main__':
   
   mjd_day = 'aug10'
   #path_file = 'transients/'+str(mjd)+'/*.txt'
   path_file = 'transients/'+mjd_day+'/*.txt'
   path_save = 'corr_transients/'+mjd_day+'/'

   try:
        os.mkdir(path_save)
   except OSError:
        print 'directry already exists'

   files_txt = sorted(glob.glob(path_file))
   for filetxt in files_txt:
       print filetxt
       inp = np.loadtxt(filetxt)

       inp_med = np.zeros(inp.shape)

       k = 0
       for i in range(inp.shape[0]):
           if search_vlss_ini(inp[i,0],inp[i,1]) == 0:
              inp_med[k,:] = inp[i,:]
              k += 1
       inp_med = inp_med[:k,:]

       out  =  filter_snr(inp_med,mjd_day)
       np.savetxt(path_save+os.path.basename(filetxt), out, '%10.3f')
