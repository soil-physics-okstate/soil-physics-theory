#PSP_diffusion1D
from __future__ import print_function, division
import numpy as np
from PSP_ThomasAlgorithm import ThomasBoundaryCondition

waterDensity = 1000.           
bulkSoilDensity = 2650.       
area = 1                       
Theta = 0.01                   
D_0   = 1e-9                  
a_const = 2.8
klQ = 0.                   
	

n = 100                        
z = np.zeros(n+2, float)       
vol = np.zeros(n+2, float)    
a = np.zeros(n+2, float)      
b = np.zeros(n+2, float)       
c = np.zeros(n+2, float)       
d = np.zeros(n+2, float)       

dz = np.zeros(n+2, float)     
Conc = np.zeros(n+2, float)   
oldConc = np.zeros(n+2, float)    
concWat = np.zeros(n+2, float)   
D = np.zeros(n+2, float)      

fd = np.zeros(n+2, float)       
fc = np.zeros(n+2, float)       
epsilon = np.zeros(n+2, float)  
        
def initializeDiffusion(n_, z_,dz_,vol_, theta):
    # vector depth [m]
    for i in range(n+2): 
        z[i] = z_[i]
        dz[i] = dz_[i]
        vol[i] = vol_[i]
    for i in range(1, n+2):
        if(z[i]>0.05 and z[i]<0.2):
           concWat[i] = 1e-3
        else:
           concWat[i] = 0
        oldConc[i] = concWat[i] * (klQ + (waterDensity*theta[i]) / bulkSoilDensity)
        Conc[i]    = oldConc[i]
            
def cellCentFiniteVolWater(dt, waterFluxDensity, theta, meanType, factor):
    g = 1.0 - factor
    fd[0] = 0
    fc[0] = 0
    for i in range(1, n+1):
        theta_mean = 1/2*(theta[i]+theta[i+1])
        D[i]   = D_0 * a_const * np.power(theta_mean,3.) \
                 + (Theta * waterFluxDensity[i]) / (waterDensity*theta_mean)
        fd[i] = area * D[i] / dz[i]         
        fc[i] = area * waterFluxDensity[i]  
        if(waterFluxDensity[i] < 0):   
              epsilon[i] = 1
    for i in range(1, n+1):    
        a[i] = -fd[i-1]*factor - fc[i-1]*(1-epsilon[i-1])*factor
        if (i == 1):
            b[i] =((klQ + (waterDensity*theta[i])/bulkSoilDensity) * vol[i]/dt 
				     + fd[i]*factor + fc[i]*(1-epsilon[i])*factor - fc[i-1]*epsilon[i-1]*factor)
            c[i] = -fd[i]*factor + fc[i]*epsilon[i]*factor
            d[i] = (oldConc[i] * vol[i]/dt + g * ((fd[i]+fc[i]*(1-epsilon[i])+
				   fc[i-1]*epsilon[i-1])*oldConc[i] - fc[i]*epsilon[i]*oldConc[i+1]))
        elif (i < n):
            b[i] =(((klQ + (waterDensity*theta[i])/bulkSoilDensity) * vol[i]/dt + 
				  fd[i-1]*factor + fd[i]*factor + fc[i]*(1-epsilon[i])*factor) 
				     - fc[i-1]*epsilon[i-1]*factor)
            c[i] = -fd[i]*factor + fc[i]*epsilon[i]*factor
            d[i] = (oldConc[i] * vol[i]/dt + g * ((fd[i-1]+fc[i-1]*(1-epsilon[i-1]))*oldConc[i-1] - 
				   (fd[i-1]+fd[i]+fc[i]*(1-epsilon[i])-fc[i-1]*epsilon[i-1])*oldConc[i] - fc[i]*epsilon[i]*oldConc[i+1]))

        elif (i == n):
            b[i] = ((klQ + (waterDensity*theta[i])/bulkSoilDensity) * vol[i]/dt 
				         + fd[n-1]*factor + fc[i]*(1-epsilon[i])*factor - fc[i-1]*epsilon[i-1]*factor)
            c[i] = 0
            d[i] = (oldConc[i] * vol[i]/dt  + g * ((fd[i-1]+fc[i-1]*(1-epsilon[i-1]))*oldConc[i-1] - 
				   (fd[i]+fc[i]*(1-epsilon[i])+fc[i-1]*epsilon[i-1])*oldConc[i]))

    ThomasBoundaryCondition(a, b, c, d, concWat, 1, n)
    return True
