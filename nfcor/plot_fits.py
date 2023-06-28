"""
Plotting the MRA fits files

"""


import numpy as np
import glob
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import scoreatpercentile as sc

files = sorted(glob.glob('/leo/savin/ovro-lwa/low_res/nfcor/mra5/sub5*.fits'))
print (files)

freq_list = ['31.296 MHz', '36.528 MHz', '41.760 MHz', '46.992 MHz']

#freq_list = ['31.296 MHz', '36.528 MHz', '46.992 MHz']
#time_list = ['13 s', '26 s', '39 s', '52 s', '65 s', '78 s']
time_list = ['13 s', '26 s', '39 s']
tm = 0

fig  = plt.figure(layout = "constrained", figsize = (14,4))
#fig  = plt.figure(layout = "constrained", figsize = (14,12))
#fig.suptitle("MRA 5")
fig.supxlabel("J2000 RA [degrees]")
fig.supylabel("J2000 Dec [degrees]")

#Ellipse locations
#MRA4 = 307,34
#MRA1 = 162,28
#MRA2 = 164,33
#MRA3 = 275, 33
#MRA5 = 60, 37

for i,filename in enumerate(files):

    hdu_list = fits.open(filename)
    hdu_list.info()
    
    #data =  hdu_list[0].data.T
    data =  hdu_list[0].data
    
    ##Normalize the data for displaying in the log space
    #vmin = data.min()
    #vmax = data.max()
    #data = (data - vmin)+1.0 #/(vmax - vmin)
    #data = 10*np.log(data)
    #vmin1 = data.min()
    #vmax1 = data.max()
     
    #data = (data - vmin1)/(vmax1 - vmin1)
    

    vmax = sc(data, 99.9)
    vmin = sc(data, 0.1)

    header = hdu_list[0].header
    print (data.shape)
    wcs_im = WCS(header)
    print (header['BMAJ'], header['BMIN'], header['BPA'])
    
    #ax = plt.subplot(1,4,i+1, projection = wcs_im)
    
    ax = plt.subplot(1,4,i+1, projection = wcs_im)

    #im = plt.imshow(data[120:350,70:180], cmap = 'jet', origin = 'lower', aspect = 'equal')
    
    im = ax.pcolormesh(data, cmap = 'jet', vmin = vmin, vmax = vmax)
            
    #im1 = ax.pcolormesh(10*np.log10(data), cmap = 'gray')
    
    #ax.axes.xaxis.set_visible(True)
    #ax.axes.yaxis.set_visible(False)
    
    ell = Ellipse((59,38), width = header['BMIN'], height = header['BMAJ'], angle = -header['BPA'], facecolor = None, edgecolor = 'white', transform = ax.get_transform('icrs'))
    ax.add_artist(ell)
    ax.tick_params(direction = 'in')
    
    ax.set_xlabel(' ')# ('J2000 RA (degrees)')
    ax.set_ylabel(' ') #('J2000 Dec (degrees)')
    
    if i < 4:
       ax.set_title(freq_list[i])
    
    if (i+1) % 4 == 1:
       ax.set_ylabel(time_list[tm])
       tm += 1
    """

    if i < 3:
       ax.set_title(freq_list[i])

    if (i+1) % 3 == 1:
       ax.set_ylabel(time_list[tm])
       tm += 1
    """
    

    cbar = plt.colorbar(im)
    cbar.set_label('Jy/beam')
    #cbar = plt.colorbar(im1)
    #cbar.set_label('[dB]')
    
    #overlay = ax.get_coords_overlay('icrs')
    #overlay.grid(color = 'white', ls = 'dotted')

#plt.figtext(0.02, 0.8, '13 s')
#plt.figtext(0.02, 0.5, '26 s')
#plt.figtext(0.02, 0.2, '39 s')
#plt.tight_layout(pad = 0.3)



#plt.savefig("mra1_grid_log.png", dpi = 150)
plt.savefig("mra5_grid.png", dpi = 150)
plt.show()    
