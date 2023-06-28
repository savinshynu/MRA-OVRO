"""
Functions used in the transient search pipeline

"""

import numpy as np
from lsl import astro

def maxi_sub(data,x_in,y_in):


    x_in=int(x_in)
    y_in=int(y_in)
    maxi=data[x_in,y_in]
    x1=x_in
    y1=y_in
    for l in range(x_in-25,x_in+26,1):
        for k in range(y_in-25,y_in+26,1):
            if data[l,k] > maxi:
               maxi=data[l,k]
               x1=l
               y1=k

    return x1,y1


def maxi_unsub(data,x_in,y_in):


    x_in=int(x_in)
    y_in=int(y_in)
    maxi=data[x_in,y_in]
    x1=x_in
    y1=y_in
    for l in range(x_in-15,x_in+16,1):
        for k in range(y_in-15,y_in+16,1):
            if data[l,k] > maxi:
               maxi=data[l,k]
               x1=l
               y1=k

    return x1,y1




def getcoord(x,y,sz,psize,MJD,h,m,s,lat,lon):

    LAT = lat*(np.pi/180.0)
    UT = h + (m + s/60.0)/60.0
    d = MJD - 40000 - 11544 + UT/24.0
    LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)
    sz = sz/2.0 
    x = ((x) - float(sz))*float(np.pi/180.0*psize)
    y = ((y) - float(sz))*float(np.pi/180.0*psize)
    ALT = np.arccos(((x**2.0)+(y**2.0))**(1.0/2.0))

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


def gethrz(x,y,sz,psize,MJD,h,m,s,lat,lon):

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
    return(az*(180.0/np.pi),ALT*(180.0/np.pi))



def trackcoord(ra,dec,MJD,h,m,s,size,psize,lat,lon):

    #print ra
    #print dec

    RA = ra
    LAT = lat*(np.pi/180.0)
    UT = h + (m + s/60.0)/60.0
    d = MJD - 40000 - 11544 + UT/24.0
    DEC = dec*(np.pi/180)
    LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)
    #print LST

    HA = np.mod(LST-RA,360)*(np.pi/180.0)
    #print HA
    ALT = np.arcsin(np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(HA))
    az = np.arccos((np.sin(DEC) - np.sin(ALT)*np.sin(LAT))/(np.cos(ALT)*np.cos(LAT)))
    B = np.sin(HA) >= 0.0
    C = np.sin(HA) < 0.0
    AZ = np.zeros(HA.shape)
    AZ[B] = 2*np.pi - az[B]
    AZ[C] = az[C]
    #print AZ

    x = np.cos(ALT)*np.sin(-AZ)
    y = np.cos(ALT)*np.cos(AZ)

    x = (x*((180.0/np.pi)/psize)) + (size/2.0)  # + 1.0)
    y = (y*((180.0/np.pi)/psize)) + (size/2.0)  # + .365)#+ 1.365)

    Att = ALT
    #print x    
    if x.shape[0] > 1:
       ind = np.where((Att < 0))
       x[ind] = 1.0
       y[ind] = 1.0

    else:
    
        if Att < 0:
           x = 1.0
           y = 1.0

    return (x,y,Att)


def mask_image(im,time,psize,lat,lon):
    time = np.reshape(time,(4,1))
    dat_im = im
    # mask for the horizon
    #print "masking"
    A = np.ones((dat_im.shape[0],dat_im.shape[1]))
    a = np.arange(dat_im.shape[0])
    b = np.arange(dat_im.shape[1])
    ag, bg = np.meshgrid(a,b,indexing='ij')
    ag = ag.astype('float')
    bg = bg.astype('float')
    #print "start calculating"
    hrz_cal  = (ag-(A.shape[0]/2.0))**2 + (bg-(A.shape[0]/2.0))**2
    #print "looking for values"
    hr_ang = 0
    hrz_lim = ((180*np.cos(hr_ang*(np.pi/180.0)))/(np.pi*psize))**2
    cond = hrz_cal > hrz_lim
    A[cond] = 0.0


    sp =  astro.get_solar_equ_coords(time[0,0] + (time[1,0] + (time[2,0] + time[3,0]/60.0)/60.0)/24.0  + 2400000.5)
    Sunra = sp[0]
    Sundec = sp[1]

    badra = np.array((83.633,350.85,299.868,187.706,139.52375,252.78375,Sunra,69.28,49.94,187.3,123.4,211.2,150.5,176.3,147,24.42,261.2,277.4,140.3,247.2))
    baddec = np.array((22.0145,58.815,40.734,12.391,-12.0956,4.9925,Sundec,29.6,41.51,2.049,48.21,51.96,28.79,31.26,7.421,33.16,-0.9939,48.75,45.67,39.55))

    for source in xrange(badra.shape[0]):

        x,y,Att = trackcoord(badra[source],baddec[source],time[0,:],time[1,:],time[2,:],time[3,:],A.shape[0],psize,lat,lon)

        if (np.isnan(x) == False) and (np.isnan(y) == False):
           x = int(np.round(x))
           y = int(np.round(y))
           xa, ya  = maxi_sub(dat_im,x,y)

           if Att > 0:
              dat_im[xa-85:xa+85,ya-85:ya+85] = 0

    dat_im_nohrz = dat_im*A
    return dat_im,dat_im_nohrz



def ang_dist_mult(x1,y1,x2,y2):
    a = np.sin(y1*(np.pi/180.0))*np.sin(y2*(np.pi/180.0))
    b = np.cos(y1*(np.pi/180.0))*np.cos(y2*(np.pi/180.0))*np.cos((x1-x2)*(np.pi/180.0))
    z=a+b
    bad = np.where((np.abs(z) > 1.0) )[0]
    if len(bad) > 0:
       for i in bad:
           if abs((z[i]-int(z[i]))/int(z[i])) < 0.01:
              z[i] = int(z[i])
               
    ang_dist = np.arccos(z)
    return ang_dist*(180.0/np.pi)


def ang_dist_sing(x1,y1,x2,y2):
    a = np.sin(y1*(np.pi/180.0))*np.sin(y2*(np.pi/180.0))
    b = np.cos(y1*(np.pi/180.0))*np.cos(y2*(np.pi/180.0))*np.cos((x1-x2)*(np.pi/180.0))
    z=a+b
    if (z >1.0) or (z<-1.0):
       if abs((z-int(z))/int(z)) < 0.01:
          z = int(z)
    ang_dist = np.arccos(z)
    return ang_dist*(180.0/np.pi)

