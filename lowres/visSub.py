"""
Subtract noise like visibilities from a peak event visibilities

"""


import sys
import aipy
import numpy

from casacore.tables import table

from astropy.constants import c as vLight
vLight = vLight.to('m/s').value

#from lsl.imaging import utils
import nf_utils

from matplotlib import pyplot as plt

# Parse the command line
#distance = float(sys.argv[1])*1e3
filename1 = sys.argv[1]
filename2 = sys.argv[2]

## Location of the MRA in latitude (deg), longitude (deg), and height (m)
#MRA = (37.80055526, -118.24152502, 100e3)

# Load the data and make a copy (we'll need it later)
ms1 = nf_utils.CorrelatedData(filename1)
ms2 = nf_utils.CorrelatedData(filename2)
data1 = ms1.get_data_set(1)
data2 = ms2.get_data_set(1)

"""
# Find the location on the sky and the position of 'distance' along 
# that line-of-sight
#mra_az, mra_el, mra_dist = ms.station.get_pointing_and_distance(MRA) # rad, rad, m
az, el, dist = [191.175141, 45.006377, 100e+3]
mra_az, mra_el, mra_dist = [az*(numpy.pi/180.0), el*(numpy.pi/180.0), dist]
mra_enz = numpy.array([numpy.sin(mra_az)*numpy.cos(mra_el),
                       numpy.cos(mra_az)*numpy.cos(mra_el),
                       numpy.sin(mra_el)])
mra_enz *= dist
mra_ra, mra_dec = ms1.station.radec_of(mra_az, mra_el)
mra_bdy = aipy.amp.RadioFixedBody(mra_ra*1.0, mra_dec*1.0)

## Re-phase to get the (u,v,w) coordinates
#data1.rephase(mra_bdy)
#data2.rephase(mra_bdy)


# Actually phase to the MRA
data1.uvw = data0.uvw*1.0
for pds in data1:
    pds.data = pds.data.astype(numpy.complex128)
    
    delays = {}
    ## Loop over baselines
    for k in range(data1.nbaseline):
        ### Load in the data
        i,j = data1.baselines[k]
        vis = pds.data[k,:]
        
        ### Get the antenna delays
        try:
           di = delays[i]
        except KeyError:
           di = numpy.sqrt(  (ms.station.antennas[i].stand.x - mra_enz[0])**2 \
                           + (ms.station.antennas[i].stand.y - mra_enz[1])**2 \
                           + (ms.station.antennas[i].stand.z - mra_enz[2])**2 )
           di = di / 0.299792	# ns
           delays[i] = di
        try:
           dj = delays[j]
        except KeyError:
           dj = numpy.sqrt(  (ms.station.antennas[j].stand.x - mra_enz[0])**2 \
                           + (ms.station.antennas[j].stand.y - mra_enz[1])**2 \
                           + (ms.station.antennas[j].stand.z - mra_enz[2])**2 )
           dj = dj * (1e9 / vLight)   # ns
           delays[j] = dj
           
        ### Apply the new phasing        
        vis = data1.antennaarray.unphs2src(vis, 'z', j, i)
        vis *= numpy.exp(2j*numpy.pi*(di-dj)*data1.antennaarray.get_afreqs())
        pds.data[k,:] = vis
"""
# Save the new information back to the measurement set
ms3 = table(filename1, readonly=False, ack=False)
print ("listing existing columns")
print (ms3.colnames())
if 'CORRECTED_DATA' in ms3.colnames():
    corname = 'CORRECTED_DATA'
else:
    corname = 'DATA'

for i,b in enumerate(data1.baselines):
    if i % 5000 == 0:
        print(i)
        
    dataT = ms3.query("ANTENNA1==%i AND ANTENNA2==%i" % b)
    uvw = dataT.getcol('UVW')
    #vis = dataT.getcol('CORRECTED_DATA')
    vis = dataT.getcol(corname)

    new_uvw = data1.uvw[i,:,0] * vLight / data1.freq[0]
    new_uvw.shape = (1,)+new_uvw.shape
    
    new_vis = numpy.array([data1.XX.data[i]-data2.XX.data[i],
                           data1.XY.data[i]-data2.XY.data[i],
                           data1.YX.data[i]-data2.YX.data[i],
                           data1.YY.data[i]-data2.YY.data[i]])
    new_vis = new_vis.T*1.0
    new_vis.shape = (1,)+new_vis.shape
    
    dataT.putcol('UVW', new_uvw)
    dataT.putcol(corname, new_vis)
ms3.close()
