#PSP_infiltration1D
from __future__ import print_function, division

import PSP_grid as grid
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
from PSP_soil import *

waterDensity = 1000.        
area = 1                    
maxNrIterations = 100
tolerance = 1e-3
n = 100                         
z = np.zeros(n+2, float)      
vol = np.zeros(n+2, float)    
a = np.zeros(n+2, float)       
b = np.zeros(n+2, float)       
c = np.zeros(n+2, float)       
d = np.zeros(n+2, float)       

dz = np.zeros(n+2, float)      
psi = np.zeros(n+2, float)     
dpsi = np.zeros(n+2, float)   
theta = np.zeros(n+2, float)   
oldTheta = np.zeros(n+2, float)  
C = np.zeros(n+2, float)      
k = np.zeros(n+2, float)      
waterFluxDensity = np.zeros(n+2, float)      

u = np.zeros(n+2, float)      
du = np.zeros(n+2, float)      
f = np.zeros(n+2, float)      
H = np.zeros(n+2, float)      
H0 = np.zeros(n+2, float)      
        
def initializeWater(funcType, soil, theta_0, solver):
    global z
    
    z = grid.linear(n, soil.lowerDepth)
    vol[0] = 0
    for i in range(n+1): 
        dz[i] = z[i+1]-z[i]
        if (i > 0): vol[i] = area * (z[i+1] - z[i-1]) / 2.0
    
    if (solver == NEWTON_RAPHSON_MFP):
        psi_0 = MFPFromTheta(soil, theta_0)
        k_0 = hydraulicConductivityFromMFP(soil, psi_0)
    else:
        psi_0 = waterPotential(funcType, soil, theta_0)
        k_0 = hydraulicConductivityFromTheta(funcType, soil, theta_0)

    psi[0] = 0
    for i in range(1, n+2):
        oldTheta[i] = theta_0
        theta[i] = theta_0
        psi[i] = psi_0
        H[i] = psi[i] - z[i]*g
        k[i] = k_0
            
def NewtonRapsonMFP(funcType, soil, dt, ubPotential, isFreeDrainage):
    b3 = (2.0 * soil.Campbell_b + 3.0) / (soil.Campbell_b + 3.0)
    #apply upper boundary condition
    airEntry = airEntryPotential(funcType, soil)
    ubPotential = min(ubPotential, airEntry) 
    psi[1] = MFPFromPsi(soil, ubPotential)
    oldTheta[1] = thetaFromMFP(soil, psi[1])
    theta[1] = oldTheta[1]
    k[1] = hydraulicConductivityFromMFP(soil, psi[1])
    psi[0] = psi[1]
    k[0] = 0.0
    
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    
    nrIterations = 0
    massBalance = 1
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        massBalance = 0
        for i in range(1, n+1):
            k[i] = hydraulicConductivityFromMFP(soil, psi[i])
            capacity = theta[i] / ((soil.Campbell_b + 3.0) * psi[i])
            C[i] = waterDensity * vol[i] * capacity / dt
            u[i] = g * k[i]
            f[i] = (psi[i+1] - psi[i]) / dz[i] - u[i]
            if (i == 1): 
                a[i] = 0
                c[i] = 0
                b[i] = 1 / dz[i] + C[i] + g * b3 * k[i] / psi[i]
                d[i] = 0
            else:
                a[i] = -1 / dz[i-1] -g * b3 * k[i-1] / psi[i-1]
                c[i] = -1 / dz[i]
                b[i] = 1 / dz[i-1] + 1 / dz[i] + C[i] + g * b3 * k[i] / psi[i]
                d[i] = f[i-1] - f[i] + (waterDensity * vol[i] 
                                        * (theta[i] - oldTheta[i]) / dt)
                massBalance += abs(d[i])
    
        ThomasBoundaryCondition(a, b, c, d, dpsi, 1, n)
        
        for i in range(1, n+1):
            psi[i] -= dpsi[i]
            psi[i] = min(psi[i], soil.CampbellMFP_he)
            theta[i] = thetaFromMFP(soil, psi[i])
            
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
            
        nrIterations += 1

    if (massBalance < tolerance):
        flux = -f[1]
        for i in range(1,n+1):
          waterFluxDensity[i] = -f[i]
        CourantNr = max(waterFluxDensity[1:np.size(waterFluxDensity)-1]/(0.5*theta[1:np.size(theta)-1]
                                            +0.5*theta[2:np.size(theta)]))*dt/min(dz[1:np.size(dz)-1])
        if(CourantNr<1.0):
           return True, nrIterations, flux
        else:
           return False, nrIterations, 0
    else:
        return False, nrIterations, 0


