import numpy as np
from matplotlib import pyplot as plt

lc1 = np.loadtxt("lc_mra1_subt.txt")
lc2 = np.loadtxt("lc_mra2_subt.txt")
lc3 = np.loadtxt("lc_mra3_subt.txt")
lc4 = np.loadtxt("lc_mra4_subt.txt")
lc5 = np.loadtxt("lc_mra5_subt.txt")


# 8,8,11,14,8 are the peak of integrations of MRA1-5
xr1 = np.arange(len(lc1))*13.0 - (3*13.0)
xr2 = np.arange(len(lc2))*13.0 - (5*13.0)
xr3 = np.arange(len(lc3))*13.0 - (10*13.0)
xr4 = np.arange(len(lc4))*13.0 - (12*13.0)
xr5 = np.arange(len(lc5))*13.0 - (7*13.0)

plt.rcParams.update({'font.size': 14})

plt.plot(xr1, lc1,  marker='o', markerfacecolor = "none",  linestyle='solid', label = "MRA1")
plt.plot(xr2, lc2,  marker='o', markerfacecolor = "none", linestyle='solid', label = "MRA2")
plt.plot(xr3, lc3,  marker='o', markerfacecolor = "none",  linestyle='solid', label = "MRA3")
plt.plot(xr4, lc4,  marker='o', markerfacecolor = "none", linestyle = 'solid',  label = "MRA4")
plt.plot(xr5, lc5, marker='o', markerfacecolor = "none", linestyle='solid', label = "MRA5")

plt.ylabel("Flux Density [Jy]")
plt.xlabel("Time [s]")
plt.ylim(-200,2400)
plt.legend()
plt.show()
