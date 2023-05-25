import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

dat1 = np.load('spec/spec_mra1.txt.npy')
dat2 = np.load('spec/spec_mra2.txt.npy')
dat3 = np.load('spec/spec_mra3.txt.npy')
dat4 = np.load('spec/spec_mra4.txt.npy')
dat5 = np.load('spec/spec_mra5.txt.npy')

freq_chan = np.linspace(dat1[0,0], dat1[-1,0], 100)

# For fitting the spectrum
def func(x,a,b):

    return (a*(x**b))

plt.rcParams.update({'font.size': 13})

plt.subplot(2,3,1)

popt, pcov = curve_fit(func,dat1[:,0],dat1[:,1],maxfev=2000)
perr = np.sqrt(np.diag(pcov))

#print popt
#print perr

plt.scatter(dat1[:,0], dat1[:,1],s = 10, marker = 'o', facecolors = 'none',edgecolors = 'b', label= 'Data')
plt.plot(freq_chan,func(freq_chan, *popt),label='SI = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Flux Density [Jy]')
plt.title("MRA1")
plt.legend()
        

plt.subplot(2,3,2)

popt, pcov = curve_fit(func,dat2[:,0],dat2[:,1],maxfev=2000)
perr = np.sqrt(np.diag(pcov))

#print popt
#print perr

plt.scatter(dat2[:,0], dat2[:,1], s = 10, marker = 'o', facecolors = 'none',edgecolors = 'b', label= 'Data')
plt.plot(freq_chan,func(freq_chan, *popt),label='SI = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Flux Density [Jy]')
plt.title("MRA2")
plt.legend()


plt.subplot(2,3,3)

popt, pcov = curve_fit(func,dat3[:,0],dat3[:,1],maxfev=2000)
perr = np.sqrt(np.diag(pcov))

#print popt
#print perr

plt.scatter(dat3[:,0], dat3[:,1], s = 10, marker = 'o', facecolors = 'none',edgecolors = 'b', label= 'Data')
plt.plot(freq_chan,func(freq_chan, *popt),label='SI = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Flux Density [Jy]')
plt.title("MRA3")
plt.legend()


plt.subplot(2,3,4)

popt, pcov = curve_fit(func,dat4[:,0],dat4[:,1],maxfev=2000)
perr = np.sqrt(np.diag(pcov))

#print popt
#print perr

plt.scatter(dat4[:,0], dat4[:,1], s = 10, marker = 'o', facecolors = 'none',edgecolors = 'b', label= 'Data')
plt.plot(freq_chan,func(freq_chan, *popt),label='SI = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Flux Density [Jy]')
plt.title("MRA4")
plt.legend()

plt.subplot(2,3,5)

popt, pcov = curve_fit(func, dat5[:,0], dat5[:,1],maxfev=2000)
perr = np.sqrt(np.diag(pcov))

#print popt
#print perr

plt.scatter(dat5[:,0], dat5[:,1], s = 10, marker = 'o', facecolors = 'none',edgecolors = 'b', label= 'Data')
plt.plot(freq_chan,func(freq_chan, *popt),label='SI = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Flux Density [Jy]')
plt.title("MRA5")
plt.legend()



plt.show()

