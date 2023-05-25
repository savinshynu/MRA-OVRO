#!/usr/bin/env python

import os
import sys
import numpy

#from casacore.tables import table

filename = sys.argv[-1]
    
print "Filename: %s" %(filename)
    
tb.open(filename, nomodify=False)
uvw = tb.getcol('UVW')
flg = tb.getcol('FLAG')
print uvw.shape
print  "Mean FLAG value: %f" % (100*flg.sum()/flg.size)
    
uvd = numpy.sqrt(uvw[0,:]**2 + uvw[1,:]**2)
wgt = 1 - flg.mean(axis=0).mean(axis=0)
print uvd.shape,wgt.shape
print  "  Mean (u,v) distance: %f " % ((uvd*wgt).sum()/wgt.sum())
tb.close()
