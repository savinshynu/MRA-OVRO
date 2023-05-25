import numpy as np
import sys
from lsl import astro
import glob
import os
#import time as tm
import h5py
from utils import mask_image, maxi_sub, ang_dist_sing, ang_dist_mult
#from scipy.stats import scoreatpercentile as sc
#from matplotlib import pyplot as plt

def getcoord(x,y,sz,psize,MJD,h,m,s,lat,lon):

    LAT = lat*(np.pi/180.0)
    UT = h + (m + s/60.0)/60.0
    d = MJD - 40000 - 11544 + UT/24.0
    LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)

    sz = sz/2.0 
    x = ((x) - float(sz))*float(np.pi/180.0*psize)
    y = ((y) - float(sz))*float(np.pi/180.0*psize)

    ALT = np.arccos(((x**2.0)+(y**2.0))**(1.0/2.0));

    x = np.reshape(np.array(x + 0.01),(-1,))
    y = np.reshape(np.array(y + 0.01),(-1,))
    az = np.zeros(x.shape)
    for i in xrange(x.shape[0]):
        if y[i] >= 0.0 and x[i] < 0.0:
            az[i] = np.arctan(y[i]/x[i]) +float(np.pi/2)
        elif y[i] < 0.0 and x[i] < 0.0:
            az[i] = np.arctan(y[i]/x[i]) + float(np.pi/2)
        elif y[i] < 0.0 and x[i] >= 0.0:
            az[i] = np.arctan(y[i]/x[i]) + float(3*np.pi/2)
        elif y[i] >= 0.0 and x[i] >= 0.0:
            az[i] = np.arctan(y[i]/x[i]) + float(3*np.pi/2)

    az = np.mod(az,2*np.pi)

    DEC = np.arcsin(np.cos(az)*(np.cos(ALT)*np.cos(LAT)) + (np.sin(ALT)*np.sin(LAT)))
    HA = np.real(np.arccos((np.sin(ALT) - (np.sin(DEC)*np.sin(LAT)))/(np.cos(DEC)*np.cos(LAT))))
    RA = np.zeros(x.shape)

    for i in xrange(x.shape[0]):
        if x[i] < 0.0:
            RA[i] = np.mod(LST[i] + HA[i]/(np.pi/180.0),360.0)
        else:
            RA[i] = np.mod(LST[i] - HA[i]/(np.pi/180.0),360.0)

    dec = DEC*(180.0/np.pi)
    ra = RA

    return(ra,dec)


def trackcoord(ra,dec,MJD,h,m,s,size,psize,lat,lon):

    RA = ra
    LAT = lat*(np.pi/180.0)
    UT = h + (m + s/60.0)/60.0
    d = MJD - 40000 - 11544 + UT/24.0
    DEC = dec*(np.pi/180)
    LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)

    HA = np.mod(LST-RA,360)*(np.pi/180.0)
    #print HA
    ALT = np.arcsin(np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(HA))
    az = np.arccos((np.sin(DEC) - np.sin(ALT)*np.sin(LAT))/(np.cos(ALT)*np.cos(LAT)))
    B = np.sin(HA) >= 0.0
    C = np.sin(HA) < 0.0
    AZ = np.zeros(HA.shape)
    AZ[B] = 2*np.pi - az[B]
    AZ[C] = az[C]

    x = np.cos(ALT)*np.sin(-AZ)
    y = np.cos(ALT)*np.cos(AZ)

    x = (x*((180.0/np.pi)/psize)) + (size/2.0)
    y = (y*((180.0/np.pi)/psize)) + (size/2.0)

    Att = ALT
    x[Att<0] = 1.0
    y[Att<0] = 1.0


    return (x,y,Att)


def cluster(ra,dec,snr):
    ra_list = ra
    dec_list = dec
    # Eliminate repeating coordinates
    ra_list_new = np.zeros(ra_list.shape)
    dec_list_new = np.zeros(dec_list.shape)
    rep_coord = []
    num_r = 0
    for j in range(ra_list.shape[0]):
        if [ra_list[j],dec_list[j]] not in rep_coord:

           ra_in = ra_list[j]*np.ones(ra_list.shape[0])
           dec_in = dec_list[j]*np.ones(ra_list.shape[0])
           diff = ang_dist_mult(ra_in,dec_in,ra_list,dec_list)
           #print diff
           mult = np.where((diff[np.logical_not(np.isnan(diff))] < 3))[0]
           #print ra_list[j],dec_list[j]
           if len(mult) == 1:
              ra_list_new[num_r] = ra_list[j]
              dec_list_new[num_r] = dec_list[j]
              num_r += 1
           elif len(mult) > 1 :
              mk = np.argmax(snr[mult]) 
              ra_list_new[num_r] = ra_list[mult[mk]]
              dec_list_new[num_r] = dec_list[mult[mk]]
          
              for l in mult:
                  rep_coord += [[ra_list[l],dec_list[l]]]
              num_r += 1
              


    ra_list_new = ra_list_new[:num_r]
    dec_list_new = dec_list_new[:num_r]
    
    # finding the peak of the cluster
    peak_ra = np.zeros(ra_list_new.shape)
    peak_dec = np.zeros(dec_list_new.shape)
    peak_snr = np.zeros(dec_list_new.shape)
    for x in range(ra_list_new.shape[0]):
        ra_ini = ra_list_new[x]*np.ones(ra.shape[0])
        dec_ini = dec_list_new[x]*np.ones(dec.shape[0])
        dist = ang_dist_mult(ra_ini,dec_ini,ra,dec)
        good = np.where((dist[np.logical_not(np.isnan(dist))]<3.0))[0]
        #print ra_list_new[x],dec_list_new[x]
        #print good
        #print snr[good]
        pk = np.argmax(snr[good])
        #print pk
        peak_ra[x] = ra[good[pk]]#np.mean(ra[good])
        peak_dec[x] = dec[good[pk]]#np.mean(dec[good])
        peak_snr[x] = snr[good[pk]]
        #print peak_snr[x]
    return peak_ra,peak_dec,peak_snr

#Searching the database with the vlss catalogue and removing sources within 3 degreees and flux greater than 50 jy
vlss = np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt")                        
def search_vlss(ra_list,dec_list,snr_list,time_list):
    ra_new = np.zeros(ra_list.shape)
    dec_new = np.zeros(dec_list.shape)
    snr_new = np.zeros(snr_list.shape)
    time_new = np.zeros(time_list.shape)
    num_new = 0
    for x in range(ra_list.shape[0]):
        ra_in = np.ones(vlss.shape[0])*ra_list[x]
        dec_in = np.ones(vlss.shape[0])*dec_list[x]

        diff = abs(ang_dist_mult(ra_in,dec_in,vlss[:,0],vlss[:,1]))
        par = np.where((vlss[:,2] > 50.0) & (diff < 3.0))[0]
        #print len(par)
        if len(par) == 0:
           ra_new[num_new] = ra_list[x]
           dec_new[num_new] = dec_list[x]
           snr_new[num_new] = snr_list[x]
           time_new[num_new,:] = time_list[x,:]
           num_new += 1
    
    if num_new > 0:   
       return ra_new[:num_new],dec_new[:num_new],snr_new[:num_new],time_new[:num_new,:]
    else:
       return np.array([]),np.array([]),np.array([]),np.array([])



def find_events(data2,time2,dur,sig_thresh,ngrid,psize,lat,lon):
    ra = []
    for i in range(data2.shape[0]):
        if i > 1:
           #print i 
           sub_dat = data2[i,:,:] - np.mean(data2[np.arange(i-dur,i,1),:,:],axis = 0)
           #print "subtracted"
           im_mask, im_mask_nohrz = mask_image(sub_dat,time2[i,:],psize,lat,lon)
           mean = np.mean(im_mask_nohrz[im_mask_nohrz != 0])
           sigma = np.std(im_mask_nohrz[im_mask_nohrz != 0])
           #print mean, sigma
           event = np.where(((im_mask_nohrz -mean) > sig_thresh*sigma) & (im_mask_nohrz != 0))
           xar = event[0]
           yar = event[1]
           #print xar
           #print yar
           #vmin = sc(im_mask_nohrz,1)
           #vmax = sc(im_mask_nohrz,99.5)
           #plt.imshow(np.transpose(im_mask_nohrz),cmap ='jet',vmin = vmin, vmax = vmax,origin ='lower')
           #plt.show()
           #print len(xar)
           if len(xar) > 0:
              time_ini = np.tile(time2[i,:],(xar.shape[0],1)) 
              snr_ini =  (im_mask_nohrz[xar,yar]-mean)/sigma 
              ra_ini, dec_ini = getcoord(xar,yar,ngrid,psize,time_ini[:,0],time_ini[:,1],time_ini[:,2],time_ini[:,3],lat,lon)
              
              b = (180.0/np.pi)*np.arcsin(np.sin(dec_ini*(np.pi/180.0))*np.cos(62.6*(np.pi/180.0)) - np.cos(dec_ini*(np.pi/180.0))*np.sin(ra_ini*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.sin(62.6*(np.pi/180.0)))
              l = 33.0 +(180.0/np.pi)*np.arcsin((np.cos(dec_ini*(np.pi/180.0))*np.sin(ra_ini*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.cos(62.6*(np.pi/180.0)) + np.sin(dec_ini*(np.pi/180.0))*np.sin(62.6*(np.pi/180.0)))/np.cos(b*(np.pi/180.0)))

              gm = np.logical_or(np.logical_and(np.logical_and(b > -10 ,b < 10), snr_ini < 8),np.logical_and(np.logical_and(b > -7 ,b < 7), snr_ini < 10))
              galmask = np.logical_not(gm)

              ra_ini = ra_ini[galmask]
              dec_ini = dec_ini[galmask]
              snr_ini = snr_ini[galmask]
              time_ini = time_ini[galmask]

              #np.savetxt('coord_info.txt',np.concatenate((np.reshape(ra_ini,(-1,1)),np.reshape(dec_ini,(-1,1)),np.reshape(snr_ini,(-1,1)),time_ini),1),'%10.3f')
              if len(ra_ini) > 0:
                 ra_cl, dec_cl, snr_cl = cluster(ra_ini,dec_ini,snr_ini)
                 time_cl = np.tile(time2[i,:],(ra_cl.shape[0],1))
                 #print ra_cl.shape
                 #print dec_cl
                 if len(ra_cl) > 0:
                    ra_ar, dec_ar, snr_ar, time_ar = search_vlss(ra_cl,dec_cl,snr_cl,time_cl)
                    #print ra_ar.shape
                    if len(ra_ar) > 0:
                       try:
                           ra = np.concatenate((ra,ra_ar),0)
                           dec = np.concatenate((dec,dec_ar),0)
                           ttime = np.concatenate((ttime,time_ar),0)
                           snr = np.concatenate((snr,snr_ar),0)

                       except NameError:
                           ra = ra_ar
                           dec = dec_ar
                           ttime = time_ar
                           snr = snr_ar

    
    if len(ra) > 0:
       return (ra,dec,ttime,snr)

    else:
       return (np.nan,np.nan,np.nan,np.nan)





# Savin, this is where we need to define what files to look at. At the moment I have it set up to just look at some file which is defined below. However when we do this on the LASI computers at LWA1 and LWA-SV the script needs to find all the files for that day and then run the script on them
something = 'False'
#path =sys.argv[1]
mjd =  58346
path = '/leo/savin/ovro-lwa/hdf5_files/aug16/*.hdf5'
files = sorted(glob.glob(path))
path_2 = '/leo/savin/ovro-lwa/transients/'+str(mjd)+'/'


for filename in files:

    print filename
    sig_thresh = 6 # sigma level for thresholding 
    Lat = 37.055166 #Station Latitude (of LWA-ovro)
    Lon = -118.28166  #Station Latitude (of LWA-ovro)

    stokes = 0 # right now we only do Stokes I, but we might want to do all 4 parameters, or at least I and V

    hf = h5py.File(filename,'r')
    dset1 = hf.get('image')
    time = hf.get('time')

    ints = dset1.shape[0] #number of integrations
    #print ints
    xSize = 4096 #size of images (x-axis)
    ySize = 4096  #size of images (x-axis)
    psize = 1.875/60.0  #angular size of a pixel (at zenith)
    
    #print "data collectd"
    for data in  [dset1]: # 
     
     ints = data.shape[0]   
     for run in xrange(3):

         if run == 0:
            runav = 1 #no averaging
            dur = 2   #subtract off the previous 2 images (26 seconds)
         elif run == 1:
              runav = 2 #average 2 integrations together (26 seconds)
              dur = 2   #subtract off previous 2 images (52 seconds)
         elif run == 2:
              runav = 5 # average up 5 integrations (65 seconds)
              dur = 1    # subtract off previous image( 65 seconds)

         nints = int(ints/runav)
        
         if nints > 3*dur: #if there are at least 3 times the number of images as the subtracted amount
 
             if runav !=1: # only do the for loop if we need to integrate
                 data2 = np.zeros((nints,data.shape[1],data.shape[2]),dtype = np.float32)
                 time2 = np.zeros((nints,time.shape[1]), dtype = np.float32)

                 for i in xrange(nints): #integrate up
                    
                     data2[i,:,:] = np.mean(data[np.arange(i*runav,((i+1)*runav )),:,:],0)
                     time2[i,:] = np.mean(time[np.arange(i*runav,((i+1)*runav)),:],0)
                
             else:
                 data2 = data
                 time2 = time
             
             ra,dec,ttime,snr = find_events(data2,time2,dur,sig_thresh,xSize,psize,Lat,Lon)

             if np.any(np.logical_not(np.isnan(ra))):
 
                 something = 'True'
                
                 if run == 0:

                     try:
                         RA5 = np.concatenate((RA5,ra),0)
                         DEC5 = np.concatenate((DEC5,dec),0)
                         TIME5 = np.concatenate((TIME5,ttime),0)
                         SNR5 = np.concatenate((SNR5,snr),0)
                     except NameError: 
                         RA5 = ra
                         DEC5 = dec
                         TIME5 = ttime
                         SNR5 = snr

                 elif run ==1:
                    
                     try:
                         RA15 = np.concatenate((RA15,ra),0)
                         DEC15 = np.concatenate((DEC15,dec),0)
                         TIME15 = np.concatenate((TIME15,ttime),0)
                         SNR15 = np.concatenate((SNR15,snr),0)
                     except NameError:
                         RA15 = ra
                         DEC15 = dec
                         TIME15 = ttime
                         SNR15 = snr

                 elif run ==2:

                     try:
                         RA60 = np.concatenate((RA60,ra),0)
                         DEC60 = np.concatenate((DEC60,dec),0)
                         TIME60 = np.concatenate((TIME60,ttime),0)
                         SNR60 = np.concatenate((SNR60,snr),0)
                     except NameError:
                         RA60 = ra
                         DEC60 = dec
                         TIME60 = ttime
                         SNR60 = snr

         else:
             print 'Not enough images in db'


if something:

    try:
        os.mkdir(path_2)
    except OSError:
        print 'directry already exists'


    try:

        np.savetxt(path_2 + str(int(mjd)) + '_13s.txt',np.concatenate((np.reshape(RA5,(-1,1)),np.reshape(DEC5,(-1,1)),TIME5,np.reshape(SNR5,(-1,1))),1),'%10.3f')

    except NameError:
        print 'no 13s transients'

    try:

        np.savetxt(path_2 + str(int(mjd))+ '_26s.txt',np.concatenate((np.reshape(RA15,(-1,1)),np.reshape(DEC15,(-1,1)),TIME15,np.reshape(SNR15,(-1,1))),1),'%10.3f')

    except NameError:

        print 'no 26s transients'


    try:

        np.savetxt(path_2 + str(int(mjd))+ '_65s.txt',np.concatenate((np.reshape(RA60,(-1,1)),np.reshape(DEC60,(-1,1)),TIME60,np.reshape(SNR60,(-1,1))),1),'%10.3f')

    except NameError:

        print 'no 65s transients'


       
