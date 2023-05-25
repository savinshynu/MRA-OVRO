import numpy as np
import matplotlib.pyplot as plt
import sys


             
#n_chan_adp = 33
tot_chan = 109

"""
def check_adp(data):
    good_adp =[]
    check = np.zeros(n_chan_adp)

    for chan in np.arange(0,tot_chan,n_chan_adp):

        if not np.array_equal(data[chan:chan+n_chan_adp],check):
           good_adp += range(chan,chan+n_chan_adp)
    return good_adp
"""

def main(data,clip=3):
    
    width = 250e+3    
    chan_width = 25e+3
    sig_old = 1
    y_dat = data
    temp_ind = np.array(range(data.shape[0]))
    ind = np.array(range(data.shape[0]))
    #ini_flag = check_adp(y_dat)
    #y_dat = y_dat[ini_flag]
    #temp_ind = temp_ind[ini_flag]
    it = 0    
    while True:
      
          spec = y_dat
          smth = spec*0.0

          # Calculate the median window size - the target is determined from the 
          # 'width' keyword, which is assumed to be in Hz
          winSize = int(1.0*width/chan_width)
          winSize += ((winSize+1)%2)
          #print  winSize
          # Compute the smoothed bandpass model
          for i in xrange(smth.size):
              mn = max([0, i-winSize/2])
              mx = min([i+winSize/2+1, smth.size])
              smth[i] = np.median(spec[mn:mx])    
  
          #plt.plot(xrange(spec.shape[0]),spec)
          #plt.plot(xrange(smth.shape[0]),smth)
          #plt.show()       
                

          diff = (spec-smth) 
          sig_mn = np.std(diff)
          sigma= sig_mn
          mean = np.mean(diff)
          opt = mean
          par = abs((sigma -sig_old)/sig_old)          
        
          if it == 0:
             good = np.argwhere(abs(diff-opt) < clip*sigma)
             temp_ind = temp_ind[good[:,0]]
             y_dat = y_dat[good[:,0]]
             sig_old = sigma
             

          elif it > 0:  
             if par > 0.3:      
                good = np.argwhere(abs(diff-opt) < clip*sigma)
                temp_ind = temp_ind[good[:,0]] 
                y_dat = y_dat[good[:,0]] 
                sig_old = sigma
                
                
             else:
                 break
             
          it += 1         
 
    #bad_chan = np.int_((np.array(list(set(ind)-set(temp_ind)))))
    #print bad_chan

    #plt.plot(temp_ind,data[temp_ind])
    #plt.show()
    
    return  temp_ind 


   

