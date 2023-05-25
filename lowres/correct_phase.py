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
filename = sys.argv[1]

## Location of the MRA in latitude (deg), longitude (deg), and height (m)
#MRA = (37.80055526, -118.24152502, 100e3)

# Load the data and make a copy (we'll need it later)
ms = nf_utils.CorrelatedData(filename)
data0 = ms.get_data_set(1)
#data1 = data0.copy()

# Find the location on the sky and the position of 'distance' along 
# that line-of-sight
#mra_az, mra_el, mra_dist = ms.station.get_pointing_and_distance(MRA) # rad, rad, m
az, el, dist = [191.175141, 45.006377, 100e+3]
mra_az, mra_el, mra_dist = [az*(numpy.pi/180.0), el*(numpy.pi/180.0), dist]
mra_enz = numpy.array([numpy.sin(mra_az)*numpy.cos(mra_el),
                       numpy.cos(mra_az)*numpy.cos(mra_el),
                       numpy.sin(mra_el)])
mra_enz *= dist
mra_ra, mra_dec = ms.station.radec_of(mra_az, mra_el)
mra_bdy = aipy.amp.RadioFixedBody(mra_ra*1.0, mra_dec*1.0)

# Re-phase to get the (u,v,w) coordinates
data0.rephase(mra_bdy)

"""
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
ms2 = table(filename, readonly=False, ack=False)
print ("listing existing columns")
print (ms2.colnames())
if 'CORRECTED_DATA' in ms2.colnames():
    corname = 'CORRECTED_DATA'
else:
    corname = 'DATA'

for i,b in enumerate(data0.baselines):
    if i % 500 == 0:
        print(i)
        
    dataT = ms2.query("ANTENNA1==%i AND ANTENNA2==%i" % b)
    uvw = dataT.getcol('UVW')
    #vis = dataT.getcol('CORRECTED_DATA')
    vis = dataT.getcol(corname)

    new_uvw = data0.uvw[i,:,0] * vLight / data0.freq[0]
    new_uvw.shape = (1,)+new_uvw.shape
    
    new_vis = numpy.array([data0.XX.data[i],
                           data0.XY.data[i],
                           data0.YX.data[i],
                           data0.YY.data[i]])
    new_vis = new_vis.T*1.0
    new_vis.shape = (1,)+new_vis.shape
    
    dataT.putcol('UVW', new_uvw)
    dataT.putcol(corname, new_vis)
ms2.close()
