import sys
#import glob
import os
#import numpy as np
import shutil

def maxi_unsub(data,x_in,y_in):


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

def eq2hrz(ra,dec,MJD,h,m,s,lat,lon):

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
    
    return (AZ*(180.0/np.pi),ALT*(180.0/np.pi))

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

    Att = ALT*(180.0/np.pi)
    #x[Att<0] = 1.0
    #y[Att<0] = 1.0
    if Att < 0:
       x = 1.0
       y = 1.0


    return (x,y,Att)

def revert(ms_name):

    rname, ext = os.path.splitext(ms_name)
    revertname = "%s_revert%s" % (rname, ext)
    shutil.copytree(ms_name, revertname)
    tb.open(revertname,nomodify=False)
    c = tb.getcol('CPARAM')
    tb.putcol('CPARAM', 1.0/c)
    tb.close()


def mjd2ut(t):

    MJD = int(t)
    h = 24.0*(t - float(MJD))
    H = int(h)
    m = 60*(h - float(H))
    M = int(m)
    s = 60*(m - float(M))
    S = int(s)

    return MJD, H, M, S

# Listing the coordinates for different sources

ra_cas = 350.866417
dec_cas = 58.811778

ra_cyg = 299.868153
dec_cyg = 40.733916

ra_vir = 187.705930
dec_vir = 12.391123

ra_tau = 83.633212
dec_tau = 22.014460


# Mention the source for tracking and change based on above info for a different one
ra_cal = [ra_cas,ra_cyg,ra_vir,ra_tau]
dec_cal = [dec_cas,dec_cyg,dec_vir,dec_tau]

def ret_flux(imfile):
    ia.open(imfile)
    B = ia.getchunk()
    im  = B[:,:,0,0]
    ia.close()
    ngrid = 4096  #size of images (x-axis)
    psize = 1.875/60.0  #angular size of a pixel (at zenith)
    elev = np.zeros(len(ra_cal))
    fjy = np.zeros(len(ra_cal))
    for i in range(len(ra_cal)):
        x, y, alt = trackcoord(ra_cal[i],dec_cal[i],mjd,h,m,s,ngrid,psize,lat_ov,lon_ov)
        #print x,y, alt
        xa,ya = maxi_unsub(im,x,y)
        #print xa, ya, alt
        fjy[i] = im[xa,ya]
        elev[i] = alt
    return fjy, elev 




filename = sys.argv[-1] # passing the initial part of the filename
#ini = os.path.basename(filename)
#ini = os.path.splitext(ini)[0]+'_'
ini_split = filename.split('_')
ini = ini_split[0]+'_'+ini_split[1]+'_'

tb.open(filename+'/OBSERVATION/')
tstart, tstop = tb.getcell('TIME_RANGE',0)
tstart /= 86400
tstop /= 86400
time = (tstart+tstop)/2.0
tb.close()

tb.open(filename+'/SPECTRAL_WINDOW/')
freq_ar = tb.getcell('CHAN_FREQ')
cfreq = np.mean(freq_ar)
cfreq = np.round(cfreq/1e+6,5)
tb.close()
#print cfreq

mjd,h,m,s = mjd2ut(time)
#print mjd,h,m,s
lat_ov = 37.055166 #Station Latitude (of LWA-ovro)
lon_ov = -118.28166  #Station Latitude (of LWA-ovro)


clean(vis=ini+'af_bp.ms', imagename = ini+'afbp_image', mode='mfs', imsize=[4096,4096],cell=['1.875arcmin'],niter=0, weighting= 'briggs', robust =0, stokes='I')

imfile = ini+'afbp_image.image'

fjy, elev = ret_flux(imfile) 

flux_dict = {"cas":fjy[0],"cyg":fjy[1],"vir":fjy[2],"tau":fjy[3] }

# changed 2000 to 1000 for more peeling

if elev[0] < 5.0 or fjy[0] < 1000:
   del flux_dict['cas']
if elev[1] < 5.0 or fjy[1] < 1000:
   del flux_dict['cyg']
if elev[2] < 5.0 or fjy[2] < 1000:
   del flux_dict['vir']
if elev[3] < 5.0 or fjy[3] < 1000:
   del flux_dict['tau']

ord_dict = sorted(flux_dict.items(), key=lambda x: x[1], reverse=True)

print ord_dict

if len(ord_dict) > 0:
  
      # Generating a model of Cyg and adding into measurement set
      for src in ord_dict:
          #print src
          os.system('rm -rf '+ini+'src.cl')
          cl.done()
          if src[0] == 'cas':
             sjy = src[1]
             cl.addcomponent(flux=sjy, dir='J2000 23h23m24.000s +58d48m54.000s', index=-0.77, spectrumtype='spectral index', freq= str(cfreq)+'MHz', label='CasA')
          if src[0] == 'cyg':
             sjy = src[1]
             cl.addcomponent(flux=sjy, dir='J2000 19h59m28.357s +40d44m02.097s', index=0.085, spectrumtype='spectral index', freq=str(cfreq)+'MHz', label='CygA')

          if src[0] == 'vir':
             sjy = src[1]
             cl.addcomponent(flux=sjy, dir='J2000 12h30m49.4233s +12d23m28.043s', index=-0.856, spectrumtype='spectral index', freq=str(cfreq)+'MHz', label='VirA')
          if src[0] == 'tau':
             sjy = src[1]
             cl.addcomponent(flux=sjy, dir='J2000 05h34m31.971s +22d00m52.06s', index=-0.299, spectrumtype='spectral index', freq=str(cfreq)+'MHz', label='TauA') 
          cl.rename(ini+'src.cl')
          cl.done()

          clearcal(ini+'af_bp.ms', addmodel=True)
          ft(ini+'af_bp.ms/', complist=ini+'src.cl', usescratch=True)


          # Phase selfcal
          os.system('rm -rf '+ini+'phase_src.cal '+ini+'phase_src_revert.cal')
          gaincal(vis=ini+'af_bp.ms/', caltable=ini+'phase_src.cal',refant='34', uvrange='>10lambda', solnorm = True, calmode='p', gaintype='G')


          # Amplitude selfcal
          os.system('rm -rf '+ini+'amp_src.cal '+ini+'amp_src_revert.cal')
          gaincal(vis=ini+'af_bp.ms/', caltable=ini+'amp_src.cal',refant='34', uvrange='>10lambda', solnorm = True, calmode='a', gaintype='G',gaintable=[ini+'phase_src.cal'])
          applycal(vis = ini+'af_bp.ms', gaintable=[ini+'phase_src.cal',ini+'amp_src.cal'], interp='linear')

          # subtracting model from the corrected data
          uvsub(ini+'af_bp.ms/')
          split(vis=ini+'af_bp.ms', outputvis=ini+'af_bp_2.ms', datacolumn='corrected')    

          # revert the calibration
          revert(ini+'phase_src.cal')
          revert(ini+'amp_src.cal')
          applycal(vis=ini+'af_bp_2.ms', gaintable=[ini+'amp_src_revert.cal',ini+'phase_src_revert.cal'], interp='linear')

    
          # split the corrected data gain
          split(vis=ini+'af_bp_2.ms', outputvis=ini+'af_bp_3.ms', datacolumn='corrected')
          os.system('rm -rf '+ini+'af_bp.ms')
          os.system('mv '+ini+'af_bp_3.ms '+ini+'af_bp.ms')
          os.system('rm -rf '+ini+'af_bp_*')
    


# Imaging the subtracted data
os.system('rm -rf '+ini+'af_sub*')
clean(vis=ini+'af_bp.ms', imagename=ini+'af_sub_mfs', mode='mfs', imsize=[4096,4096],cell=['1.875arcmin'],niter=0, weighting= 'briggs', robust =0, stokes='I')

#os.system('rm -rf '+ini+'af_bp*')
#os.system('rm -rf '+ini+'afbp*')
os.system('rm -rf '+ini+'*phase*.cal')
os.system('rm -rf '+ini+'*amp*.cal')
os.system('rm -rf '+ini+'src.cl')
os.system('rm -rf '+ini+'*.flux')
os.system('rm -rf '+ini+'*.model')
os.system('rm -rf '+ini+'*.psf')
os.system('rm -rf '+ini+'*.residual')
