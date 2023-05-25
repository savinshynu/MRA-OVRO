import numpy as np
import matplotlib.pyplot as plt
from sliding_rfi_flagger_man import main
from lsl import astro
import glob
import sys
import os
from utils import mask_image,maxi_sub,maxi_unsub,trackcoord
import h5py
#from lc_correct_man import snr_in_lc
from scipy.stats import chisquare,scoreatpercentile as sc


def flag_narrow_rfi(capt):
       
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    mjd = float(sys.argv[3])
    h = float(sys.argv[4])
    m = float(sys.argv[5])
    s = float(sys.argv[6])
  
    t_in = h + m/60.0 + s/3600.00
    
    if capt not in [13,26,65]:
       sys.exit("Wrong integration length")

    if capt == 13:
       diff_tm = 6
    elif capt == 26:
       diff_tm = 20 #13
    elif capt == 65:
       diff_tm = 35

    stokes = 0

    #find the corresponding .hdf5 file for the time and .oims file fo image
    #path1 = '/leo/savin/ovro-lwa/low_res/hdf5_files/'+day+'/ovro_'+day+'_01band_'+str(format(int(h),'02d'))+'h.hdf5'
    
    path1 = sys.argv[8]
    files1  = sorted(glob.glob(path1))
    for filename1 in files1:
        hf=h5py.File(filename1,'r')
        dset1=hf.get('image')
        dset = np.mean(dset1,axis=3)
        tm_ind = hf.get('time')
        header = hf.get('header')
        ints = dset.shape[0] #number of integrations
        ngrid  = header[0,1] #size of images (x-axis)
        psize = header[0,0] #angular size of a pixel (at zenith)
        
        t_comp = tm_ind[:,1] + tm_ind[:,2]/60.0 + tm_ind[:,3]/3600.00
        try:
           peak = np.where((abs(t_comp-t_in)*3600.00 < diff_tm))[0]
        except IndexError:
           continue
        if len(peak) == 0:
           continue
        peak_ind = int(np.median(peak))
        noise = np.array(range(max(0,peak_ind-6),max(0,peak_ind-3))+range(min(ints,peak_ind+3),min(ints,peak_ind+6)))
        print peak
        
        sub_dat = np.mean(dset[peak,:,:],axis=(0)) - np.mean(dset[noise,:,:],axis=(0))
        dat_img = sub_dat
        print np.std(dat_img[1400:1800,600:900])
        time = np.array([mjd,h,m,s])

        lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
        lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)

        xp,yp,alt = trackcoord(np.array([ra]),np.array([dec]),time[0],time[1],time[2],time[3],ngrid,psize,lat_ovro,lon_ovro)
        
        xa,ya = maxi_sub(dat_img,xp,yp)
        print xp,yp,(alt*180.0/np.pi),xa,ya
        
        dat_im, dat_im_nohrz = mask_image(dat_img.copy(),time,psize,lat_ovro,lon_ovro)


        fill = np.where((dat_im_nohrz != 0))
       
        plt.imshow(np.transpose(dat_im_nohrz),cmap='jet',origin = 'lower')
        plt.colorbar()
        plt.show()   
        
        print np.std(dat_im_nohrz[fill])
        print 'IM SNR : %f,' % ((dat_im_nohrz[xa,ya]/np.std(dat_im_nohrz[fill]))) 
        
       
        hf.close()
               
flag_narrow_rfi(int(sys.argv[7]))


