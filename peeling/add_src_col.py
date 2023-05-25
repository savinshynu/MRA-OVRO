from casacore.tables import table, tableutil
import numpy as np
import sys

dir = sys.argv[-1]

ta = table(dir+'/OBSERVATION/')
tstart, tstop = ta.getcell('TIME_RANGE',0) 
ta.close()

td = table(dir+'/FIELD/')
dir_rad = td.getcell('REFERENCE_DIR',0)[0]
field_name = td.getcell('NAME',0)
td.close()


tstart /= 86400
tstop /= 86400

col1  = tableutil.makearrcoldesc('DIRECTION', 0.0, 1, comment='Direction (e.g. RA, DEC).', 
                                         keywords={'QuantumUnits':['rad','rad'], 
                                                   'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                   })
col2  = tableutil.makearrcoldesc('PROPER_MOTION', 0.0, 1, 
                                         comment='Proper motion', 
                                         keywords={'QuantumUnits':['rad/s',]})
col3  = tableutil.makescacoldesc('CALIBRATION_GROUP', 0, 
                                         comment='Number of grouping for calibration purpose.')
col4  = tableutil.makescacoldesc('CODE', "none", 
                                         comment='Special characteristics of source, e.g. Bandpass calibrator')
col5  = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         comment='Interval of time for which this set of parameters is accurate', 
                                         keywords={'QuantumUnits':['s',]})
col6  = tableutil.makescacoldesc('NAME', "none", 
                                         comment='Name of source as given during observations')
col7  = tableutil.makescacoldesc('NUM_LINES', 0, 
                                         comment='Number of spectral lines')
col8  = tableutil.makescacoldesc('SOURCE_ID', 0, 
                                         comment='Source id')
col9  = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', -1, 
                                         comment='ID for this spectral window setup')
col10 = tableutil.makescacoldesc('TIME', 0.0,
                                         comment='Midpoint of time for which this set of parameters is accurate.', 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
col11 = tableutil.makearrcoldesc('TRANSITION', 'none', 1, 
                                         comment='Line Transition name')
col12 = tableutil.makearrcoldesc('REST_FREQUENCY', 1.0, 1, 
                                         comment='Line rest frequency', 
                                         keywords={'QuantumUnits':['Hz',], 
                                                   'MEASINFO':{'type':'frequency', 
                                                               'Ref':'LSRK'}
                                                   })
col13 = tableutil.makearrcoldesc('SYSVEL', 1.0, 1, 
                                         comment='Systemic velocity at reference', 
                                         keywords={'QuantumUnits':['m/s',], 
                                                   'MEASINFO':{'type':'radialvelocity', 
                                                               'Ref':'LSRK'}
                                                   })
        
desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13])

nSource = 1
tb = table("%s/SOURCE" % dir, desc, nrow=nSource, ack=False)
        

time = (tstart+tstop)/2*86400

tb.putcell('DIRECTION', 0, dir_rad)
tb.putcell('PROPER_MOTION', 0, [0.0, 0.0])
tb.putcell('CALIBRATION_GROUP', 0, 0)
tb.putcell('CODE', 0, 'none')
tb.putcell('INTERVAL', 0, 0.0)
tb.putcell('NAME', 0, field_name)
tb.putcell('NUM_LINES', 0, 0)
tb.putcell('SOURCE_ID', 0, 0)
tb.putcell('SPECTRAL_WINDOW_ID', 0, -1)
tb.putcell('TIME', 0, time)
#tb.putcell('TRANSITION', i, [])
#tb.putcell('REST_FREQUENCY', i, [])
#tb.putcell('SYSVEL', i, [])
            
tb.close()


tc = table("%s" % dir, readonly=False, ack=False)
tname = "SOURCE"
stc = table("%s/%s" % (dir,tname), ack=False)
tc.putkeyword(tname,stc)
stc.close()
tc.close()

