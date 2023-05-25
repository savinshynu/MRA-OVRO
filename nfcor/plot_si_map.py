import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS

filename = sys.argv[1]

si_map =  np.load(filename)

print si_map.shape

# Build WCS object

wcs_dict = {
    'CTYPE1': 'RA---SIN',
    'CUNIT1': 'deg',
    'CDELT1': -0.03125,
    'CRPIX1': 1,
    'CRVAL1': 60.06685,
    'NAXIS1': 251,
    'CTYPE2': 'DEC--SIN',
    'CUNIT2': 'deg',
    'CDELT2': 0.03125,
    'CRPIX2': 1,
    'CRVAL2': 36.841193,
    'NAXIS2': 451
}
wcs_map= WCS(wcs_dict)

print wcs_map

fig = plt.figure()
ax = plt.subplot(projection=wcs_map)
plt.imshow(si_map[120:350,70:180], origin='lower', cmap='cividis', aspect='equal')
plt.xlabel(r'J2000 RA (degrees)')
plt.ylabel(r'J2000 Dec (degrees)')
cbar = plt.colorbar()
cbar.set_label('SI')
overlay = ax.get_coords_overlay('icrs')
overlay.grid(color='white', ls='dotted')

plt.show()
