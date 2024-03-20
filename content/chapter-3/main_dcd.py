# main_dcd

from math import exp
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
import PSP_grid as grid
import matplotlib.pyplot as plt
import numpy as np

def gasSolver(boundaryLayerCond, boundaryConc_O2, dg_O2, respRate, totalDepth, n):
    a  = np.zeros(n+2, float)  
    b  = np.zeros(n+2, float)  
    c  = np.zeros(n+2, float)  
    d  = np.zeros(n+2, float) 
    g  = np.zeros(n+2, float) 
    u  = np.zeros(n+2, float)  
    co = np.zeros(n+2, float)  
    
    g[0] = boundaryLayerCond
    co[0] = boundaryConc_O2
    # vector depth [m]
    z = grid.linear(n, totalDepth)
    
    # initialize matrix
    for i in range(1, n+1):
        u[i] = respRate * exp(-z[i] / 0.3) * (z[i + 1] - z[i - 1]) / 2.0
        if i < n:
            g[i] = dg_O2 / (z[i + 1] - z[i])
        else:
            g[i] = 0
        a[i + 1] = -g[i]
        b[i] = g[i - 1] + g[i]
        c[i] = -g[i]
        d[i] = u[i]

    d[1] = d[1] + g[0] * co[0]
    
    ThomasBoundaryCondition(a, b, c, d, co, 1, n)
    
    return(z, co)


def main():
    R = 8.3143                     
    n = 20                        
    totalDepth = 0.5               
    bulkDensity = 1090.            #changed from 1200 to 1090 kg/m3           
    particleDensity = 2650.         
    waterContent = 0.31            #changed from 0.2 to 0.31             
    respRate = -0.001               
    oxygenDiff = 1.77e-5           #changed from 1.39e-5 to 1.77e-5 for O2 
    temperature = 25.             
    atmPressure = 101.3           
    boundaryLayerCond = 0.001      #changed from 0.01 to 0.001
    volumeFrac_O2 = 0.21           #defined
    bg = 0.9                       #constant value ranging from 0.5-1.0 and depends on the value chosen for mg (can be seen in table 3.2.)
    mg = 2.3                       #constant value ranging from 1-2 depending on the shape of the soil particles (can be seen in table 3.2.)   
    
    porosity = 1. - bulkDensity / particleDensity
    gasPorosity = porosity - waterContent


 # O2 concentration in air [g/m^3]
    boundaryConc_O2 = (volumeFrac_O2 * atmPressure * 1000. * 32. /            #equation 3.15
                          (R * (temperature + 273.15))) 
    
 # O2 binary diffusion coefficient [m2/s]
    binaryDiffCoeff_O2 = (oxygenDiff * (101.3 / atmPressure)                  #equation 3.13
                * ((temperature + 273.15) / 273.15)**1.75)
    
    dg_O2 = binaryDiffCoeff_O2 * bg * gasPorosity**mg                         #equation 3.14
    
    z, co = gasSolver(boundaryLayerCond, boundaryConc_O2, 
                      dg_O2, respRate, totalDepth, n)
    
    po = co * R * (temperature + 273.15) / (atmPressure * 1000. * 32.)        #partial pressure of O2
                           
    pco2 = 0.21 - po                                                          #partial pressure of CO2 assuming recipirocal relationship
      
    print ("node   depth [m]   Po2  Pco2")
    for i in range(n + 2):
        print ("%3d    %6.2f      %.3f       %.3f" %(i, z[i], po[i], pco2[i]))

    print ("gas-filled porosity (cm3/cm3)")
    print("%.2f" %(gasPorosity))
    
    # plot results
    fig = plt.figure(figsize=(10,8))
    for i in range(n+1):
        plt.plot(po[i], -z[i], 'rx', pco2[i], -z[i], 'go')
    plt.legend(['O2', 'CO2'], loc="upper center", frameon = False)  
    plt.xlabel('Concentration [g m$^{-3}$]',fontsize=20,labelpad=8)
    plt.ylabel('Depth [m]',fontsize=20,labelpad=8)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=8)
    plt.show()
main()