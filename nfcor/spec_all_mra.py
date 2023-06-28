"""
Plot the spectral index maps of MRAs

"""


import numpy as np
import numpy.ma as ma
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import ndimage

keys = ['conv_1a*done*sub*', 'conv_1b*done*sub*', 'conv_1c*done*sub*']
#, 'conv_2d*done*sub*', 'conv_2e*done*sub*', 'conv_2f*done*sub*']
#keys = ['conv_4a*done*sub*', 'conv_4b*done*sub*', 'conv_4c*done*sub*']
#keys = ['conv_3a*done*sub*']

times = np.arange(len(keys))*13 
#fig, axs  = plt.subplots(2,len(keys), constrained_layout=True, figsize = (4,8))
fig, axs  = plt.subplots(2,len(keys), constrained_layout=True, figsize = (12,8))
for k,key in enumerate(keys):

    #path of image files
    path = "/leo/savin/ovro-lwa/low_res/nfcor/mra1/" + key

    #Getting image files
    files = sorted(glob.glob(path))
    print(files)

    #opening the first image to get the shape
    ia.open(files[0])
    im_shape = ia.shape()
    print(im_shape)
    ia.close()

    big_im = np.zeros((len(files), im_shape[0], im_shape[1]))
    print(big_im.shape)
    freqs = np.zeros(len(files))
    #S = np.zeros((im_shape[0], im_shape[1]))

    for i,filename in enumerate(files):
        ia.open(filename)
        im = ia.getchunk()
        print(im.shape)
        freq =  ia.summary()['refval'][-1]/1e+6
        freqs[i] = freq
        big_im[i,:,:] = im[:,:,0,0]
        ia.close()

    #Lets do the fitting part

    def fit(x, a, b):
        return a*(x**b)

    big_im = np.transpose(big_im, axes = (0,2,1))
    big_im_avg = np.mean(big_im, axis = 0)
    S = big_im_avg*0.0

    sigma = np.std(big_im_avg)
    mean = np.mean(big_im_avg)

    good = np.argwhere((np.abs(big_im_avg-mean) > 3.5*sigma))
    bad =  np.argwhere((np.abs(big_im_avg-mean) <= 3.5*sigma))



    #Masking the zero regions
    big_im_avg = ma.MaskedArray(big_im_avg)
    S = ma.MaskedArray(S)
    #big_im_avg[bad[:,0], bad[:,1]] = ma.masked
    S[bad[:,0], bad[:,1]] = ma.masked


    for j, coord in enumerate(good):
        x = coord[0]
        y = coord[1]
        dat = big_im[:,x,y]
        if j%1000 == 0:
            print(j)

        try:
            popt, pcov = curve_fit(fit, freqs, dat, maxfev = 2000)
            #print(popt)
            #perr = np.sqrt(np.diag(pcov))
            S[x,y] = popt[1]
        except RuntimeError:
            #continue
            #print (popt[1])
            S[x,y] = ma.masked

    print(type(big_im_avg))

    if len(keys) == 1:
        im1 = axs[0].pcolormesh(big_im_avg, cmap = 'hot')
        axs[0].set_title(f"{times[k]} s")
        axs[0].tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        plt.colorbar(im1, ax = axs[0], label = 'Jy/beam')

        im2 = axs[1].pcolormesh(S, cmap = 'hot')
        #axs[1].set_title("SI")
        axs[1].tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        plt.colorbar(im2, ax = axs[1], label = 'SI')

    else:
        im1 = axs[0,k].pcolormesh(big_im_avg, cmap = 'hot')
        axs[0,k].set_title(f"{times[k]} s")
        #axs[0,k].set_xticks([])
        #axs[0,k].set_yticks([])
        axs[0,k].tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        plt.colorbar(im1, ax = axs[0,k], label = 'Jy/beam')

        im2 = axs[1,k].pcolormesh(S, cmap = 'hot')
        #axs[1,k].set_title("SI")
        axs[1,k].tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        plt.colorbar(im2, ax = axs[1,k], label = 'SI')

#fig.suptitle('Average Intensity and Spectral Index')
plt.savefig("si_map_mra1")
plt.close()





