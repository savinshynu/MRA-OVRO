"""
Pipeline to conduct the transient search with OVRO low resolution image

"""

import numpy as np
import sys
from lsl import astro
import glob
import os
#import time as tm
import h5py


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




def findtranspim(data,time,dur,sigma,sz,size,runav,psize,lat,lon):
   
    tmap = np.zeros(data.shape)

    A = np.ones((data.shape[0],data.shape[1])) # mask for the horizon
    for a in xrange(A.shape[0]):
        for b in xrange(A.shape[1]):
            if (float(a)-(A.shape[0]/2.0))**2 + (float(b)-(A.shape[0]/2.0))**2 >  ((180*np.cos(30.0*(np.pi/180.0)))/(np.pi*psize))**2: #masking below 30 degrees (180* cos(theta)/(np.pi * psize)
            #if ((float(a)-(A.shape[0]/2.0))/(140.0/np.pi))**2 + ((float(b)-(A.shape[0]/2.0))/(140.0/np.pi))**2 > 1.0:
               A[a,b] = 0

    #print A.shape


    sp =  astro.get_solar_equ_coords(time[0,0] + (time[1,0] + (time[2,0] + time[3,0]/60.0)/60.0)/24.0  + 2400000.5)
    Sunra = sp[0]
    Sundec = sp[1]

    badra = np.array((83.633,350.85,299.868,187.706,139.52375,252.78375,Sunra,69.28,49.94,187.3,123.4,211.2,150.5,176.3,147,24.42,261.2,277.4,140.3,247.2))
    baddec = np.array((22.0145,58.815,40.734,12.391,-12.0956,4.9925,Sundec,29.6,41.51,2.049,48.21,51.96,28.79,31.26,7.421,33.16,-0.9939,48.75,45.67,39.55))

    #print np.concatenate((np.reshape(badra,(-1,1)),np.reshape(baddec,(-1,1))),1)


    for source in xrange(badra.shape[0]):
        
        x,y,Att = trackcoord(badra[source],baddec[source],time[0,:],time[1,:],time[2,:],time[3,:],A.shape[0],psize,lat,lon)
        #print x.shape
        #print y.shape


       #if source ==1:
        #    print x

        if source == 0:
            X = np.reshape(np.round(x),(-1,1))
            Y = np.reshape(np.round(y),(-1,1))
        else:
            X = np.concatenate((X,np.reshape(np.round(x),(-1,1))),1)
            Y = np.concatenate((Y,np.reshape(np.round(y),(-1,1))),1)
            

    
    
    #print X[:,1]
    X = np.delete(X,np.arange(0,dur),0)
    Y = np.delete(Y,np.arange(0,dur),0)
    #print X.shape
    c = np.concatenate((np.reshape(X,(X.shape[0],1,-1)),np.reshape(Y,(Y.shape[0],1,-1))),1)
    #print c[1,:,:]
    #print ken
    A = np.tile(np.reshape(A,(A.shape[0],A.shape[1],1)),(1,1,data.shape[2]))
    #print A.shape
    
    
    if runav == 0:
        tmap = ((data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)]) > np.tile(np.reshape(sigma*np.std(np.reshape(data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)],(-1,data.shape[2]-dur)),0) - np.mean((data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)]),(0,1)),(1,1,-1)),(data.shape[0],data.shape[1],1)))

    else:

        for i in xrange(dur):

#            print i 
#            print np.arange(i,data.shape[2] - (dur-i))
            if i == 0:
                datasub = data[:,:,np.arange(i,data.shape[2] - (dur-i))]
            else:
                datasub =  datasub + data[:,:,np.arange(i,data.shape[2] - (dur-i))]

    
        datasub = datasub/dur
        datasub2 = np.swapaxes(data[:,:,np.arange(dur,data.shape[2])] - datasub,0,1)
        s = np.std(np.reshape(datasub2,(-1,datasub2.shape[2])),0)
        sm = np.median(s)
        #print sm
        #print ken

#        print datasub2.shape
        

        for i in xrange(c.shape[0]):
            rc = np.argwhere(datasub2[:,:,i] > (50*np.std(datasub2[:,:,i] + np.mean(datasub2[:,:,i]))))
            rowtest = rc[:,0]
            coltest = rc[:,1]
            test = np.zeros((rowtest.shape[0]))

#            if i == 100 or i == 90:
                #print s[i]
                #print 2*sm

            
            for j in xrange(rowtest.shape[0]):
     

           
                if np.any((rowtest[j] >= (c[i,1,:]-3)) & (rowtest[j] <= (c[i,1,:]+3)) & (coltest[j] >= (c[i,0,:]-3)) & (coltest[j] <= (c[i,0,:]+3))):
                    
                    test[j] = 1
                #    print test[j]

            if np.any(test):
               # print 'yes'
                tmap[:,:,i] = 0
            elif s[i] > 2*sm:
                
                tmap[:,:,i] = 0
            else:
                for j in xrange(c.shape[2]):

                    if c[i,0,j] > 1.0:
                        
                        if (c[i,1,j]-size) < 1.0 or (c[i,0,j]-size) < 1.0 or c[i,0,j] + size > A.shape[0] or (c[i,1,j] + size) > A.shape[0]:
                            size = size - 2


                        maskx = np.arange(int(np.round(c[i,1,j]-size)),int(np.round(c[i,1,j] + size + 1.0)))
                        masky = np.arange(int(np.round(c[i,0,j]-size)),int(np.round(c[i,0,j] + size + 1.0)))
                        
                        for x in maskx:
                            for y in masky:
                                datasub2[x,y,i] = 0.0


                D = datasub2[:,:,i]
                

                S = np.std(D[D!=0.0])
                M = np.mean(D[D!=0.0])
                D = D*A[:,:,i]


                tmap[:,:,i] = D > np.tile(np.reshape(sigma*S + M,(1,1)),(data.shape[0],data.shape[1]))


                if np.any((A[:,:,i] == 0) & (datasub2[:,:,i]> np.tile(np.reshape(2*sigma*S + M,(1,1,-1)),(data.shape[0],data.shape[1])))):
                    tmap[:,:,i] = 0

                if np.std(data[:,:,i])==0 and np.std(data[:,:,i+dur]- data[:,:,i]) ==0:
                    tmap[:,:,i] = 0

                            
                
    tt = np.argwhere(np.sum(tmap,(0,1)) > 0.0)[:,0]

 #   print tt
    tmap = tmap[:,:,tt]
    tdata = data[:,:,tt + dur]
    tmdata = data[:,:,tt]

#    
    ttime = time[:,tt+dur]
    #print ttime.shape

    
    for k in xrange(tdata.shape[2]):
        
        rc = np.argwhere(tmap[:,:,k] > 0.0)

        
        row = rc[:,0]
        col = rc[:,1]


        CC = (tdata[:,:,k] - tmdata[:,:,k])*A[:,:,0]

        rc = np.argwhere(CC == np.max(CC))
        row2 = rc[:,0]
        col2 = rc[:,1]


        corr = np.round(np.sqrt((row.astype('float')+1)**2.0 + (col.astype('float')+1)**2.0)/6.0)
        ucorr = np.unique(np.round(np.sqrt((row.astype('float')+1)**2.0 + (col.astype('float')+1)**2.0)/6.0))


        signifk = np.zeros((ucorr.shape[0]))
        urow = np.zeros((ucorr.shape[0]))
        ucol = np.zeros((ucorr.shape[0]))


        for rd in xrange(ucorr.shape[0]):
            signifk[rd] = np.max(CC[row[corr == ucorr[rd]],col[corr == ucorr[rd]]])/np.std(CC)
            urow[rd] = np.mean(row[corr ==ucorr[rd]])
            ucol[rd] = np.mean(col[corr ==ucorr[rd]])

            

        corr2 = np.round(np.sqrt((row2.astype('float') + 1.0)**2.0 + (col2.astype('float') + 1.0)**2.0)/6.0)
        ucorr2 = np.unique(np.round(np.sqrt((row2.astype('float')+1.0)**2.0 + (col2.astype('float') + 1.0)**2.0)/6.0))

        signifk2 = np.zeros((ucorr2.shape[0]))
        urow2 = np.zeros((ucorr2.shape[0]))
        ucol2 = np.zeros((ucorr2.shape[0]))

        for rd in xrange(ucorr2.shape[0]):
            signifk2[rd] = np.max(CC[row2[corr2 == ucorr2[rd]],col2[corr2 == ucorr2[rd]]])/np.std(CC)
            urow2[rd] = np.mean(row2[corr2 ==ucorr2[rd]])
            ucol2[rd] = np.mean(col2[corr2 ==ucorr2[rd]])


        if np.std(row) + np.std(col) < 20:
          #  print psize
    #        print ucol
    #        print urow
            rak,deck = getcoord(ucol,urow,sz,psize,np.tile(ttime[0,k],ucol.shape),np.tile(ttime[1,k],ucol.shape),np.tile(ttime[2,k],ucol.shape),np.tile(ttime[3,k],ucol.shape),lat,lon)

            if np.any(np.isnan(rak)):
                print urow[np.isnan(rak)]
                print ucol[np.isnan(rak)]

        
        else:
            rak,deck = getcoord(ucol2,urow2,sz,psize,np.tile(ttime[0,k],ucol2.shape),np.tile(ttime[1,k],ucol2.shape),np.tile(ttime[2,k],ucol2.shape),np.tile(ttime[3,k],ucol2.shape),lat,lon)

            signifk = signifk2
            if np.any(np.isnan(rak)):
                print urow2[np.isnan(rak)]
                print ucol2[np.isnan(rak)]





        if k == 0:

            ra = rak
            dec = deck
            signif = signifk

            ttimek = np.tile(np.reshape(ttime[:,k],(4,1)),(1,np.max(rak.shape)))
            #ttimek = ttime[:,k]
            #print ttimek.shape

        else:
            
            ra = np.append(ra,rak)
            dec = np.append(dec,deck)
            signif = np.append(signif,signifk)

            ttimek = np.concatenate((ttimek,np.tile(np.reshape(ttime[:,k],(4,1)),(1,np.max(rak.shape)))),1)
     #       print ttimek.shape

            



    try:
        ttime = ttimek

    #print dec
    #print ra
        b = (180.0/np.pi)*np.arcsin(np.sin(dec*(np.pi/180.0))*np.cos(62.6*(np.pi/180.0)) - np.cos(dec*(np.pi/180.0))*np.sin(ra*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.sin(62.6*(np.pi/180.0)))
    #print np.sin(dec*(np.pi/180.0))*np.cos(62.6*(np.pi/180.0)) - np.cos(dec*(np.pi/180.0))*np.sin(ra*(np.pi/180.0) - (282.5*(np.pi/180.0))*np.sin(62.6*(np.pi/180.0))
    #print b



        l = 33.0 +(180.0/np.pi)*np.arcsin((np.cos(dec*(np.pi/180.0))*np.sin(ra*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.cos(62.6*(np.pi/180.0)) + np.sin(dec*(np.pi/180.0))*np.sin(62.6*(np.pi/180.0)))/np.cos(b*(np.pi/180.0)))

    #print l
    

        #print b
        #print l
    #print ken

#        gm = np.logical_and(np.logical_and(np.logical_and(np.logical_and(b > -10 ,b < 10), l > -50 ), l < 50), signif < 8)

        gm = np.logical_or(np.logical_and(np.logical_and(b > -10 ,b < 10), signif < 8),np.logical_and(np.logical_and(b > -7 ,b < 7), signif < 10))



   #     print gm

        galmask = np.logical_not(gm)
   # print galmask.shape
   # print ttime.shape
  #      print ra
  #      print dec


        ra = ra[galmask]
        dec = dec[galmask]
        ttime = ttime[:,galmask]
        
        return (ra,dec,ttime)
    
    except UnboundLocalError:

        return (np.nan,np.nan,np.nan)



# Savin, this is where we need to define what files to look at. At the moment I have it set up to just look at some file which is defined below. However when we do this on the LASI computers at LWA1 and LWA-SV the script needs to find all the files for that day and then run the script on them
something = 'False'

mjd = 'aug16' 
path = 'avg_hdf5/'+mjd+'/*h.hdf5'
files = sorted(glob.glob(path))
path_2 ='transients/'+mjd+'/'



for filename in files:

    print filename
    sig = 6 # sigma level for thresholding 
    source_size = 10 #estimate of the number of pixels that need to be subtracted off surrounding a source
    
    Lat = 37.055166 #Station Latitude (of LWA-ovro)
    Lon = -118.28166  #Station Latitude (of LWA-ovro)



    hf=h5py.File(filename,'r')
    dset1 = hf.get('image')
    dset2 = hf.get('time')
    dset3 = hf.get('header')
    hdr  =  dset3[0,:]
    ints = dset1.shape[0] #number of integrations
    xSize = int(hdr[4]) #size of images (x-axis)
    ySize = int(hdr[4])  #size of images (x-axis)
    psize = hdr[3] #angular size of a pixel (at zenith)

    #time = np.zeros((4,ints),dtype = np.float32)
    #data = np.zeros((xSize,ySize,ints),dtype=np.float32)
    #stokes = 0 # right now we only do Stokes I, but we might want to do all 4 parameters, or at least I and V



    #for i in xrange(ints):


        #data[:,:,i] = dset1[i,stokes,:,:]
        #time[:,i] = np.reshape(dset2[i,:],(4,))


    data = np.transpose(dset1,axes = (1,2,0))
    time = np.transpose(dset2)
    
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
                data2 = np.zeros((data.shape[0],data.shape[1],nints))
                time2 = np.zeros((time.shape[0],nints))

                for i in xrange(nints): #integrate up
                    
                    data2[:,:,i] = np.mean(data[:,:,np.arange(i*runav,((i+1)*runav ))],2)
                    time2[:,i] = np.mean(time[:,np.arange(i*runav,((i+1)*runav))],1)
                
            else:
                data2 = data
                time2 = time

            ra,dec,ttime = findtranspim(data2,time2,dur,sig,xSize,source_size,runav,psize,Lat,Lon)

            if np.any(np.logical_not(np.isnan(ra))):

                something = 'True'
                
                if run == 0:

                    try:
                        RA5 = np.concatenate((RA5,ra),0)
                        DEC5 = np.concatenate((DEC5,dec),0)
                        TIME5 = np.concatenate((TIME5,ttime),1)
                      #  print RA5.shape
                      #  print TIME5.shape

                    except NameError: 
                        RA5 = ra
                        DEC5 = dec
                        TIME5 = ttime
                elif run ==1:
                    
                    try:
                        RA15 = np.concatenate((RA15,ra),0)
                        DEC15 = np.concatenate((DEC15,dec),0)
                        TIME15 = np.concatenate((TIME15,ttime),1)
                    except NameError:
                        RA15 = ra
                        DEC15 = dec
                        TIME15 = ttime
                elif run ==2:

                    try:
                        RA60 = np.concatenate((RA60,ra),0)
                        DEC60 = np.concatenate((DEC60,dec),0)
                        TIME60 = np.concatenate((TIME60,ttime),1)
                    except NameError:
                        RA60 = ra
                        DEC60 = dec
                        TIME60 = ttime


      
        else:
            print 'Not enough images in db'


if something:

    try:
        os.mkdir(path_2)
    except OSError:
        print 'directry already exists'


    try:

        np.savetxt(path_2 + mjd + '_13s.txt',np.concatenate((np.reshape(RA5,(1,-1)).T,np.reshape(DEC5,(1,-1)).T,TIME5.T),1),'%10.2f')

    except NameError:
        print 'no 13s transients'

    try:

        np.savetxt(path_2 + mjd+ '_26s.txt',np.concatenate((np.reshape(RA15,(1,-1)).T,np.reshape(DEC15,(1,-1)).T,TIME15.T),1),'%10.2f')

    except NameError:

        print 'no 26s transients'


    try:

        np.savetxt(path_2 + mjd+ '_65s.txt',np.concatenate((np.reshape(RA60,(1,-1)).T,np.reshape(DEC60,(1,-1)).T,TIME60.T),1),'%10.2f')

    except NameError:

        print 'no 65s transients'
