import numpy as np
import numpy.ma as ma
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import ndimage

#path of image files
path = "/leo/savin/ovro-lwa/low_res/nfcor/mra3/conv_3a*done*sub*"

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

good = np.argwhere((np.abs(big_im_avg-mean) > 3.0*sigma))
bad =  np.argwhere((np.abs(big_im_avg-mean) <= 3.0*sigma))


#plt.pcolormesh(big_im_avg.T, cmap = "jet")
#plt.savefig("big_im2")
#plt.show()
#plt.close()


#Masking the zero regions
big_im_avg = ma.MaskedArray(big_im_avg)
S = ma.MaskedArray(S)
big_im_avg[bad[:,0], bad[:,1]] = ma.masked
S[bad[:,0], bad[:,1]] = ma.masked


for i, coord in enumerate(good):
    x = coord[0]
    y = coord[1]
    dat = big_im[:,x,y]
    if i%1000 == 0:
        print(i)
        #print((i/good.shape[0])*100)

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

#Enter the angle in which to rotate
ang = -45.0

big_im_avg_rot = ndimage.rotate(big_im_avg, ang, reshape = False)
S_rot = ndimage.rotate(S, ang, reshape = False)
S_rot_avg = np.mean(S, axis = 1)


fig, (ax1, ax2, ax3)  = plt.subplots(1,3, constrained_layout=True, figsize = (12,8))

im1 = ax1.pcolormesh(big_im_avg, cmap = 'hot')
ax1.set_title("Averaged map")
plt.colorbar(im1, ax = ax1)

im2 = ax2.pcolormesh(S, cmap = 'hot')
ax2.set_title("Spectral Index map")
plt.colorbar(im2, ax = ax2)

ax3.plot(S_rot_avg, '.')
ax3.set_title("Spectral Index")
ax3.set_ylabel("SI")
ax3.set_xlabel("Pixels")

plt.savefig("si_map_3a")
plt.close()


fig, (ax1, ax2, ax3)  = plt.subplots(1,3, constrained_layout=True, figsize = (12,8))

im1 = ax1.pcolormesh(big_im_avg_rot, cmap = 'hot')
ax1.set_title("Averaged map")
plt.colorbar(im1, ax = ax1)

im2 = ax2.pcolormesh(S_rot, cmap = 'hot')
ax2.set_title("Spectral Index map")
plt.colorbar(im2, ax = ax2)

ax3.plot(S_rot_avg, '.')
ax3.set_title("Spectral Index")
ax3.set_ylabel("SI")
ax3.set_xlabel("Pixels")

plt.savefig("si_map_rot_3a")
plt.close()



