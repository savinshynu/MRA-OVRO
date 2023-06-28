"""
Clustering sources from candidate list
"""

import numpy as np
import matplotlib.pyplot as plt
from utils import ang_dist_sing, ang_dist_mult
from scipy.cluster.vq import kmeans2

data = np.loadtxt("coord_info.txt")
ra = data[:,0]
dec = data[:,1]

#plt.scatter(ra,dec)
#plt.show()


ra_list = np.array(ra[0])
dec_list = np.array(dec[0])

num_k = 0
old_ra = ra[0]
old_dec =  dec[0]

for i in range(ra.shape[0]):
    if i>0:
       cur_ra = ra[i]
       cur_dec = dec[i]
       sep = ang_dist_sing(cur_ra,cur_dec,old_ra,old_dec)
       if sep > 3:
          ra_list = np.append(ra_list,cur_ra)
          dec_list = np.append(dec_list,cur_dec)
          old_ra = cur_ra
          old_dec = cur_dec
          num_k += 1

print num_k

ra_list_new = np.zeros(ra_list.shape)
dec_list_new = np.zeros(dec_list.shape)
rep_coord = []
#rep_dec = np.array([])
num_r = 0
for j in range(ra_list.shape[0]):
    ra_in = ra_list[j]*np.ones(ra_list.shape[0])
    dec_in = dec_list[j]*np.ones(ra_list.shape[0])
    diff = ang_dist_mult(ra_in,dec_in,ra_list,dec_list)
    mult = np.where((diff[np.logical_not(np.isnan(diff))] < 3))[0]
    #print mult
    if len(mult) == 1:
       ra_list_new[num_r] = ra_list[j]
       dec_list_new[num_r] = dec_list[j]
       num_r += 1
    elif len(mult) > 1 :
       if [ra_list[j],dec_list[j]] not in rep_coord:
          ra_list_new[num_r] = np.mean(ra_list[mult])
          dec_list_new[num_r] = np.mean(dec_list[mult])
          #print ra_list[j],dec_list[j]
          #print rep_coord
          for l in mult: 
              rep_coord += [[ra_list[l],dec_list[l]]]
          num_r += 1
ra_list_new = ra_list_new[:num_r]
dec_list_new = dec_list_new[:num_r]

print num_r
#for x in range(ra_list.shape[0]):
#    print ra_list[x], dec_list[x]

print "mean values"
mean_ra = np.zeros(ra_list_new.shape)
mean_dec = np.zeros(dec_list_new.shape)
for x in range(ra_list_new.shape[0]):
    ra_ini = ra_list_new[x]*np.ones(ra.shape[0])
    dec_ini = dec_list_new[x]*np.ones(dec.shape[0])
    dist = ang_dist_mult(ra_ini,dec_ini,ra,dec)
    good = np.where((dist[np.logical_not(np.isnan(dist))]<3.0))[0]
    #good = np.where((dist[np.logical_not(np.isnan(dist))]<3.0))[0]
    ra_mean = np.mean(ra[good])
    dec_mean = np.mean(dec[good])
    print ra_list_new[x],dec_list_new[x] #,ra_mean,dec_mean


#ini = np.concatenate((np.reshape(ra_list_new,(-1,1)),np.reshape(dec_list_new,(-1,1))),1)
#centroid, label = kmeans2(data[:,:2], ini, minit='matrix')

#print centroid

plt.scatter(ra,dec)
plt.show()

