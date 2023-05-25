import sys
import os,glob
import numpy as np
import matplotlib.pyplot as plt

f = open('/users/savin/ovro_meteor/utctimes.txt','r')
fh  = f.readlines()
list1 = []
day  = []
time = []
for line in fh:
    day_line = line.split(" ")[0]
    time_line = line.split(" ")[1]
    filename = line.split(" ")[3][:43]
    list1.append(filename)
    day.append(day_line)
    time.append(time_line)
    
    #print line
f.close()

"""
for h in range(18,24):
    path  = '/wheeler/scratch/savin/%02d/*.ms' % (h)
    files = sorted(glob.glob(path))

    f2 = open('/wheeler/scratch/savin/time_10aug_%02dh.txt' % (h),'w')

    for i,filename in enumerate(files):
        file_base  = os.path.basename(filename)
        file_ext  = os.path.splitext(file_base)[0]
        ind =  list1.index(file_ext)
        f2.write("%d  %s  %s \n" % (i,day[ind],time[ind]))
"""

def trackcoord(ra,dec,MJD,h,m,s,lat,lon):

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

ha = []
ma = []
sa = []
name = []

for i,x in enumerate(day):

    if x == '2018-08-16':
       tm = time[i]
       hi = int(tm.split(':')[0])
       mi = int(tm.split(':')[1])
       si = int(tm.split(':')[2])
       ha.append(hi)
       ma.append(mi)
       sa.append(si)
       name.append(list1[i])

ha = np.array(ha)
ma = np.array(ma)
sa = np.array(sa)
mjd = 58346

ra_cas = 350.866417
dec_cas = 58.811778

ra_cyg = 299.868153
dec_cyg = 40.733916

ra_vir = 187.705930
dec_vir = 12.391123

ra_tau = 83.633212
dec_tau = 22.014460

#ra_ar = np.ones(len(ha))*ra_cyg
#dec_ar = np.ones(len(ha))*dec_cyg
#mjd_ar = np.ones(len(ha))*mjd

lat_ov = 37.055166 #Station Latitude (of LWA-ovro)
lon_ov = -118.28166  #Station Latitude (of LWA-ovro)

az, alt = trackcoord(ra_cyg, dec_cyg, mjd, ha, ma, sa, lat_ov, lon_ov )

plt.plot(range(alt.shape[0]),alt)
plt.show()

print name[1718],ha[1718]
