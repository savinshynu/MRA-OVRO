import sys
import glob
import os
import numpy as np
import shutil
import time as tm
import h5py

def mjd2ut(t):

    MJD = int(t)
    h = 24.0*(t - float(MJD))
    H = int(h)
    m = 60*(h - float(H))
    M = int(m)
    s = 60*(m - float(M))
    S = int(s)

    return MJD, H, M, S


path1_base = sys.argv[-1]
path1 = path1_base+'/*.ms'
files = sorted(glob.glob(path1))
#hour = path1_base.split('/')[-1]
#day = path1_base.split('/')[-2]
day = 'aug16'
hour = '13'
sub = '07'

#path2 = '/users/savin/wheeler-scratch/'+day
path2 = path1_base[:-3]

# bandpass calibration file already derived 
path3 = '/users/savin/peeling_ovro/aug16/'
cal_file = path3+sub+'_calibrationtable.bcal'
print cal_file
 
capt = len(files)
if capt == 0:
   sys.exit('No files in the directory')

os.system('rm -rf '+path2+'/ovro_%s_%sband_%sh.hdf5'%(day,sub,hour))
outfile = h5py.File(path2+'/ovro_%s_%sband_%sh.hdf5'%(day,sub,hour),'a')
# Initializing a hdf5  dataset
dset = outfile.create_dataset('image',(capt,256,256,109),dtype='float32') #just stokes I
tset = outfile.create_dataset('time',(capt,4),dtype='float32') #just 
hset = outfile.create_dataset('header',(capt,5),dtype='float32') #just

for intg,filename in enumerate(files):

    #t0 = tm.time()
    print filename
        
    ini = os.path.splitext(filename)[0]+'_'
    #print ini
    
    try:
       tb.open(filename+'/OBSERVATION/')
       tstart, tstop = tb.getcell('TIME_RANGE',0)
       tstart /= 86400
       tstop /= 86400
       time = (tstart+tstop)/2.0
       tb.close()
    except RuntimeError:
       continue

    mjd,h,m,s = mjd2ut(time)
    #print mjd,h,m,s
    
    try:
       tb.open(filename+'/SPECTRAL_WINDOW/')
       freq_ar = tb.getcell('CHAN_FREQ')
       #cfreq = np.mean(freq_ar)
       #cfreq = np.round(cfreq/1e+6,5)
       freq_ini = freq_ar[0]/1e+6
       freq_last = freq_ar[-1]/1e+6
       freq_width = (freq_ar[1]-freq_ar[0])/1e+6
       tb.close()
    except RuntimeError:
       continue

    # Flagging bad channels, antennas and baselines    
    #flagdata(vis = filename, mode='manual', spw='0:000;076')
    flagdata(vis = filename, mode='manual', antenna = '9,77,82,87,89,90,91,92,93,104,107,114,115,116,117,118,119,120,121,127,128,145,148,153,157,161,164,168,179,197,201,203,219,220,222,224,225,236,238,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255')
    flagdata(vis = filename, mode='manual', antenna ='0&1;0&2;1&2;1&3;2&3;2&4;3&4;3&5;4&5;4&6;5&6;5&7;6&7;6&8;6&58;7&8;8&10;10&11;10&12;11&12;11&13;12&13;12&14;13&14;13&15;14&15;14&16;15&16;15&17;15&58;16&17;16&18;17&18;17&19;18&19;18&20;19&20;19&21;20&21;20&22;21&22;21&23;22&23;22&24;22&26;23&24;23&25;24&25;24&26;25&26;25&27;26&27;26&28;27&28;27&29;28&29;28&30;29&30;29&31;30&31;30&32;31&32;31&33;32&33;32&34;33&34;33&35;34&35;34&36;35&36;35&37;35&68;36&37;36&38;37&38;37&39;38&39;38&40;39&40;39&41;40&41;40&42;41&42;41&43;42&43;42&44;43&44;43&45;43&46;43&62;43&123;44&45;44&46;45&46;45&47;46&47;46&48;47&48;47&49;48&49;48&50;48&51;49&50;49&51;50&51;50&52;51&52;51&53;52&53;52&54;53&54;53&55;54&55;54&56;55&56;55&57;55&85;56&57;56&58;57&58;57&59;57&191;58&59;58&60;58&63;58&135;59&60;59&61;60&61;60&62;61&62;61&63;62&63;62&64;62&80;62&123;62&149;62&188;62&223;62&233;63&64;63&65;63&123;63&135;63&188;64&65;64&66;65&66;65&67;66&67;66&68;67&68;67&69;68&69;68&70;69&70;69&71;69&100;69&103;70&71;70&72;71&72;71&73;71&106;72&73;72&74;73&74;73&75;74&75;74&76;75&76;76&78;78&79;78&80;79&80;79&81;80&81;80&122;80&190;81&83;81&85;83&84;83&85;83&86;84&85;84&86;85&86;86&88;88&122;88&190;94&95;94&96;95&96;95&97;96&97;96&98;97&98;97&99;98&99;98&100;99&100;99&101;100&101;100&102;100&103;101&102;101&103;102&103;103&105;105&106;106&108;108&109;108&110;109&110;109&111;109&144;110&111;110&112;111&112;111&113;112&113;122&123;122&124;122&190;123&124;123&125;123&188;123&191;123&199;124&125;124&126;125&126;125&130;125&191;129&130;129&131;130&131;130&132;131&132;131&133;132&133;132&134;133&134;133&135;134&135;134&136;134&163;135&136;135&137;136&137;136&138;136&165;137&138;137&139;138&139;138&140;139&140;139&141;140&141;140&142;141&142;141&143;142&143;142&144;143&144;144&146;146&147;146&173;147&149;149&150;149&151;150&151;150&152;150&177;150&180;151&152;151&180;152&154;154&155;154&156;155&156;156&158;158&159;158&160;159&160;160&162;160&187;162&163;163&165;165&166;165&167;166&167;167&169;169&170;169&171;169&206;170&171;170&172;171&172;171&173;172&173;172&174;173&174;173&175;174&175;174&176;175&176;175&177;176&177;176&178;177&178;178&180;180&181;180&182;181&182;181&183;182&183;182&184;182&217;183&184;183&185;184&185;184&186;184&191;185&186;185&187;185&192;185&230;185&231;186&187;186&188;187&188;187&189;187&192;187&199;187&202;187&213;187&231;188&189;188&190;188&191;189&190;189&191;190&191;190&192;191&192;191&193;191&239;192&193;192&194;193&194;193&195;194&195;194&196;195&196;196&198;198&199;198&200;198&226;199&200;200&202;202&204;204&205;204&206;204&232;205&206;205&207;206&207;206&208;206&210;207&208;207&209;208&209;208&210;209&210;209&211;210&211;210&212;211&212;211&213;211&214;212&213;212&214;213&214;213&215;214&215;214&216;215&216;215&217;216&217;216&218;217&218;221&223;226&227;226&228;227&228;227&229;228&229;228&230;229&230;229&231;230&231;230&232;231&232;231&233;232&233;232&234;233&234;233&235;234&235;235&237;237&239')

    #Applying bandpass first    
    applycal(vis = filename, gaintable=[cal_file], interp='linear')

    # Splitting out the corrected data
    os.system('rm -rf '+ini+'af_bp.ms  ' + ini+'channel*')
    split(vis = filename, outputvis = ini+'af_bp.ms', datacolumn='corrected')
   
    #t1 = tm.time()
    #print 'preprocessing time %f' %(t1-t0)
    #clean(vis=ini+'af_bp.ms', imagename = ini+'afbp_image', mode='mfs', imsize=[256,256],cell=['31.68arcmin'],uvrange = '<20lambda', niter=0, weighting= 'natural', stokes='I')
    
    clean(vis = ini+'af_bp.ms', imagename = ini+'channel', mode='channel', imsize=[256,256],cell=['30arcmin'],uvrange = '<30lambda', niter=0, weighting= 'natural', stokes='I')

    os.system('rm -rf '+ini+'af_bp*')
    os.system('rm -rf '+ini+'*.flux')
    os.system('rm -rf '+ini+'*.model')
    os.system('rm -rf '+ini+'*.psf')
    os.system('rm -rf '+ini+'*.residual')
   
    #t2 = tm.time()
    #print 'clean+rm time %f' %(t2-t1)
    tset[intg,:] = mjd,h,m,s
    hset[intg,:] = freq_ini,freq_last,freq_width,31.68/60.0,256
    #Once the .image is created, then use ia function in casa to openthat file and copy the data
    ia.open(ini+'channel.image')
    B = ia.getchunk()
    dset[intg,:,:,:] = B[:,:,0,:]
    ia.close()
    #t3 = tm.time()
    #print 'transfer time %f' %(t3-t2)
    #print intg
    continue


# closing the write out file
outfile.close()
