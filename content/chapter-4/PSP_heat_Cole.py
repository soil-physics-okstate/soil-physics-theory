#PSP_heat_Cole
from __future__ import print_function, division
import PSP_grid as grid
from PSP_ThomasAlgorithm import *
import numpy as np
from PSP_heatSoil import *

area = 1                    	
maxNrIterations = 100
tolerance = 1.e-2               

n = 20                        
z = np.zeros(n+2, float)  
zCenter = np.zeros(n+2, float)      
dz = np.zeros(n+2, float)       
vol = np.zeros(n+2, float)     
wc = np.zeros(n+2, float)       
a = np.zeros(n+2, float)        
b = np.zeros(n+2, float)        
c = np.zeros(n+2, float)        
d = np.zeros(n+2, float)       

T = np.zeros(n+2, float)      
dT = np.zeros(n+2, float)       
oldT = np.zeros(n+2, float)     
C_T = np.zeros(n+2, float)     
lambda_ = np.zeros(n+2, float)  
k_mean = np.zeros(n+2, float)   
f = np.zeros(n+2, float)       


def initialize(T_0, thetaIni, solver):
    global z, dz, zCenter, vol, wc, T, oldT
    # vector depth [m]
    z = grid.geometric(n, 1.0)
  
    vol[0] = 0
    for i in range(n+1): 
        dz[i] = z[i+1]-z[i]
        if (i > 0): vol[i] = area * dz[i]
    for i in range(n+2): 
        zCenter[i] = z[i] + dz[i]*0.5
        
    if (solver == CELL_CENT_FIN_VOL):
        for i in range(n+1): 
            dz[i] = zCenter[i+1]-zCenter[i]
                 
    for i in range(1, n+2):
        T[i] = T_0
        oldT[i] = T_0
        wc[i] = thetaIni
    return z
    

def finiteDifference(airT, boundaryT, dt, factor):
    g = 1.0 - factor
    energyBalance = 1.
    for i in range(1, n+2):
        T[i] = oldT[i]
    nrIterations = 0
    while ((energyBalance > tolerance) and (nrIterations < maxNrIterations)):
        for i in range(1, n+2):
            #lambda_[0]= 20.
            lambda_[i] = thermalConductivity(bulkDensity, wc[i], clay, T[i])
            C_T[i] = heatCapacity(bulkDensity, wc[i])*vol[i]
        f[0] = 0.
        for i in range(1, n+1):
            f[i]=area* lambda_[i] / dz[i]
        for i in range(1, n+1):    
            if (i == 1):
                a[i] = 0.
                b[i] = 1.
                c[i] = 0.
                d[i] = airT
            elif (i < n):
                a[i] = -f[i-1]*factor
                b[i] = C_T[i]/dt + f[i-1]*factor + f[i]*factor
                c[i] = -f[i]*factor
                d[i] = C_T[i]/dt * oldT[i] +(1.-factor)*(f[i-1]*oldT[i-1]+f[i]
                                            *oldT[i+1]-(f[i-1]+f[i])*oldT[i])
            elif (i == n):
                a[n] = 0.
                b[n] = 1.
                c[n] = 0.
                d[n] = boundaryT 
        ThomasBoundaryCondition(a, b, c, d, T, 1, n)
        dSum = 0
        for i in range(2, n):
            dSum += C_T[i]*(T[i]-oldT[i])
        energyBalance = (abs(dSum - factor*dt*(f[1]*(T[1]-T[2]) 
                    - f[n-1]*(T[n-1]-boundaryT)) - g*dt*(f[1]*(oldT[1]-oldT[2]) 
                    - f[n-1]*(oldT[n-1]-boundaryT))))
        nrIterations += 1
        
    if (energyBalance < tolerance):
        flux = f[1]*(T[1]-T[2])
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0

        
def cellCentFiniteVol(airT, boundaryT, dt):
    energyBalance = 1.
    for i in range(1, n+2):
        T[i] = oldT[i] 
    nrIterations = 0
    while ((energyBalance > tolerance) and (nrIterations < maxNrIterations)):
        for i in range(1, n+2):
            lambda_[i] = thermalConductivity(bulkDensity, wc[i], clay, T[i])
            C_T[i] = heatCapacity(bulkDensity, wc[i])*vol[i]
        f[0] = 0.
        for i in range(1, n+1):
            f[i] = area * kMean(LOGARITHMIC, lambda_[i], lambda_[i+1]) / dz[i]
        for i in range(1, n+1):    
            if (i == 1):
                a[i] = 0.
                b[i] = 1.
                c[i] = 0.
                d[i] = airT
            elif (i < n):
                a[i] = -f[i-1]
                b[i] = C_T[i]/dt + f[i-1] + f[i]
                c[i] = -f[i]
                d[i] = C_T[i]/dt * oldT[i]
            elif (i == n):
                a[n] = 0.
                b[n] = 1.
                c[n] = 0.
                d[n] = boundaryT 
        ThomasBoundaryCondition(a, b, c, d, T, 1, n)
        dSum = 0
        for i in range(2, n):
            dSum += C_T[i]*(T[i]-oldT[i])
        energyBalance = (abs(dSum - f[1]*(T[1]-T[2])*dt 
                           + f[n-1]*(T[n-1]-boundaryT)*dt))
        nrIterations += 1
        
    if (energyBalance < tolerance):
        flux = f[1]*(T[1]-T[2])
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0
    

