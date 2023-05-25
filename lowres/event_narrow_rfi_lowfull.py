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
from scipy.optimize import curve_fit

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
       diff_tm = 20#13
    elif capt == 65:
       diff_tm = 35

    stokes = 0

    #find the corresponding .hdf5 file for the time and .oims file fo image
    #path1 = '/leo/savin/ovro-lwa/low_res/hdf5_files/'+day+'/ovro_'+day+'_01band_'+str(format(int(h),'02d'))+'h.hdf5'
    
    path1 = sys.argv[8]
    files1  = sorted(glob.glob(path1))
    for filename1 in files1:
        print filename1
        hf=h5py.File(filename1,'r')
        dset=hf.get('image')
        tm_ind = hf.get('time')
        header = hf.get('header')
        chan = hf.get('channel') 
        #head  =  header[0,:]
        ints = dset.shape[0] #number of integrations
        ngrid  = header[0,1] #size of images (x-axis)
        psize = header[0,0] #angular size of a pixel (at zenith)
        #nchan = 109
        #start_freq, stop_freq, bandwidth = head[:3]
        #print 'CF: %f' %((start_freq+stop_freq)/2.0)
       
        #print ngrid,psize
        #print tm_ind[:,0]
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
        print noise 
        print ints
        sub_dat = np.mean(dset[peak,:,:,:],axis=(0)) - np.mean(dset[noise,:,:,:],axis=(0))
        dat_img = np.mean(sub_dat,axis=2)

        time = np.array([mjd,h,m,s])

        lat_ovro = 37.055166 #Station Latitude (of LWA-ovro)
        lon_ovro = -118.28166  #Station Latitude (of LWA-ovro)

        xp,yp,alt = trackcoord(np.array([ra]),np.array([dec]),time[0],time[1],time[2],time[3],ngrid,psize,lat_ovro,lon_ovro)
        
        xai,yai = maxi_sub(dat_img,xp,yp)
        print xp,yp,(alt*180.0/np.pi),xai,yai
        
      
        plt.plot((np.arange(sub_dat.shape[2])),sub_dat[xai,yai,:],color='blue', marker='.', linestyle='solid')
        #plt.xlabel("Frequency (MHz)")
        #plt.ylabel("Flux_response")
        plt.show()
       

        gd = main(sub_dat[xai,yai,:])
        #gd = list(set(gd)-set(range(218,248)))
        #gd = np.arange(31,73)
        dat_img_flag = np.mean(sub_dat[:,:,gd],axis=2)
        #print np.std(dat_img_flag[80:140,40:90])
        xa, ya = maxi_unsub(dat_img_flag,xai,yai)
        #xa,ya = [206,210]   
           
        src_spec = sub_dat[xa,ya,gd]

        freq = np.zeros(436)
        freq[:] = chan[0,:]  #np.linspace(start_freq,stop_freq,nchan)
        freq = freq[gd]
        #np.save('spec/spec_mra5.txt', np.concatenate((np.reshape(freq,(-1,1)), np.reshape(src_spec,(-1,1))),axis = 1))
        plt.scatter(gd,src_spec,s=20,color='blue', label= 'Data')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Flux Density (Jy)')
        #plt.plot(freq,spec_func(freq, *popt),label='Spectral index = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
        plt.show()
        
        def func(x,a,b):

            return (a*(x**b))


        popt, pcov = curve_fit(func,freq,src_spec,maxfev=3000)
        perr = np.sqrt(np.diag(pcov))
        
        #xdat = np.linspace(freq[0],freq[-1],100)
        plt.scatter(freq,src_spec,s =20,color='blue', label= 'Data')
        plt.plot(freq,func(freq, *popt),label='Spectral index = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
    
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Flux Density (Jy)')
        plt.legend()
        plt.show()
        
        
        dat_im, dat_im_nohrz = mask_image(dat_img.copy(),time,psize,lat_ovro,lon_ovro)
        dat_im_flag, dat_im_flag_nohrz = mask_image(dat_img_flag.copy(),time,psize,lat_ovro,lon_ovro)

        #mean_flag = np.mean(dat_img_flag)
        #mean_noflag = np.mean(dat_img)

        fill_flag = np.where((dat_im_flag_nohrz != 0))
        fill = np.where((dat_im_nohrz != 0))
       
        plt.imshow(np.transpose(dat_img),cmap='jet',origin='lower')
        plt.show()   
       
        plt.pcolormesh(np.transpose(dat_img_flag),cmap='jet')
        #plt.pcolormesh(np.transpose(dat_im_flag_nohrz),cmap='jet')
        plt.colorbar()
        plt.show()
        print np.std(dat_im_nohrz[fill]),np.std(dat_im_flag_nohrz[fill_flag])
        #print dat_im_nohrz[xa,ya]
        print 'IM SNR bf flag: %f, af flag: %f' % ((dat_im_nohrz[xa,ya]/np.std(dat_im_nohrz[fill])), (dat_im_flag_nohrz[xa,ya]/np.std(dat_im_flag_nohrz[fill_flag]))) 
        
        """
        if sys.argv[9] == 'c':
           ## results from the light curve
           lc_sig_noflag,lc_sig,lc_lin_pol,lc_circ_pol = snr_in_lc(ra,dec,db,peak,stokes)
           print 'LC SNR bflag: %f, SNR aflag: %f, Lpol:%f, Cpol: %f' % (lc_sig_noflag,lc_sig,lc_lin_pol,lc_circ_pol)
           #lc_sig,lc_lin_pol,lc_circ_pol = snr_in_lc(ra,dec,db,peak,stokes)
           #print lc_sig,lc_lin_pol,lc_circ_pol
        """   
       
        hf.close()
               
flag_narrow_rfi(int(sys.argv[7]))


