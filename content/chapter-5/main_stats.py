#PSP_waterRetentionFitting
from __future__ import print_function, division
try: input = raw_input
except: pass

import numpy as np
import matplotlib.pyplot as plt
from PSP_readDataFile import readDataFile
from PSP_Marquardt import *

def main_stats(file_str, model_int):
    bold = '\033[1m'
    # read experimental values
    myOutput, isFileOk = readDataFile(file_str, 1, '\t', False)
    if (not isFileOk): 
        print('Wrong file: error reading row nr.', myOutput)
        return(False)
    waterPotential = myOutput[:,0]
    waterContent = myOutput[:,1]
    
    # select water retention curve
    
    waterRetentionCurve = 0
    while (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > CAMPBELL_IPPISCH_VG):
        waterRetentionCurve = float(model_int)
        if (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > CAMPBELL_IPPISCH_VG):
            print('wrong choice.')
  
    # initialize parameters
    thetaS = max(waterContent)
    #thetaR = min(waterContent)
    thetaR = 0.08
    air_entry = 1.0
    Campbell_b = 4.0
    VG_alpha = 1/air_entry
    VG_n = 1.5
    VG_m = 1. - 1./VG_n
      
    if (waterRetentionCurve == CAMPBELL):
        b0 = np.array([thetaS, air_entry, Campbell_b], float)
        bmin = np.array([thetaS, 0.1, 0.1], float)
        bmax = np.array([thetaS*1.1, 20., 10.], float)
        print(bold + 'CAMPBELL' + bold)
    elif (waterRetentionCurve == VAN_GENUCHTEN):
        b0 = np.array([thetaS, thetaR, VG_alpha, VG_n, VG_m], float)
        bmin = np.array([thetaS, 0.0, 0.01, 0.01, 0.01], float)
        bmax = np.array([1.0, thetaR, 10., 10., 1.], float)
        print(bold + 'VAN_GENUCHTEN' + bold)
    elif (waterRetentionCurve == RESTRICTED_VG):
        b0 = np.array([thetaS, thetaR, VG_alpha, VG_n], float)
        bmin = np.array([thetaS, 0.0, 0.01, 1.], float)
        bmax = np.array([1, thetaR, 10., 10.], float)
        print(bold + 'RESTRICTED_VG' + bold)
    elif (waterRetentionCurve == IPPISCH_VG):
        b0 = np.array([thetaS, thetaR, air_entry, VG_alpha, VG_n], float)
        bmin = np.array([thetaS, 0.0, 0.1, 0.01, 1.], float)
        bmax = np.array([1, thetaR, 10., 10., 10.], float)
    elif (waterRetentionCurve == CAMPBELL_IPPISCH_VG):
        b0 = np.array([thetaS, thetaR, air_entry, VG_alpha, VG_n], float)
        bmin = np.array([thetaS, 0.0, 0.1, 0.01, 1.], float)
        bmax = np.array([1, thetaR, 10., 10., 10.], float)
   
    else:
        print ('wrong choice.')
        return(False)

    print ("\nFitting")
    b = Marquardt(waterRetentionCurve, b0, bmin, bmax, waterPotential, waterContent)

    print ("\nthetaS = ", b[0])
    if (waterRetentionCurve == CAMPBELL):
        print ("AirEntry = ", b[1])
        print ("b = ", b[2])
    elif (waterRetentionCurve == VAN_GENUCHTEN):
        print ("thetaR = ", b[1])
        print ("alpha = ", b[2])
        print ("n = ", b[3])
        print ("m = ", b[4])
    elif (waterRetentionCurve == RESTRICTED_VG):
        print ("thetaR = ", b[1])
        print ("alpha = ", b[2])
        print ("n = ", b[3])
    elif (waterRetentionCurve == IPPISCH_VG):
        print ("thetaR = ", b[1])
        print ("AirEntry = ", b[2])
        print ("alpha = ", b[3])
        print ("n = ", b[4])
    elif (waterRetentionCurve == CAMPBELL_IPPISCH_VG):
        print ("thetaR = ", b[1])
        print ("AirEntry = ", b[2])
        print ("alpha = ", b[3])
        print ("n = ", b[4])

    myWP = np.logspace(-5, 8, 500)
    myWC = estimate(waterRetentionCurve, b, myWP)


#main()    
