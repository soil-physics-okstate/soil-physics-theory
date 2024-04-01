#PSP_heatSoil.py
from __future__ import division
from math import exp
import numpy as np

GEOMETRIC = 0
LOGARITHMIC = 1

FIN_DIFF = 1
CELL_CENT_FIN_VOL = 2

# bulkDensity = 1500. # kg/m^3
# clay = 0.10 # g/g

def kMean(meanType, k1, k2):
    if (meanType == GEOMETRIC): 
        k = np.sqrt(k1 * k2)
    elif (meanType == LOGARITHMIC):
        if (k1 == k2):
            k = k1
        else:
            k = (k1-k2) / np.log(k1/k2)  
    return k

def thermalConductivity(bulkDensity, waterContent, clay, temperature):
    ga = 0.088                
    thermalConductivitysolid = 2.5  
    atmPressure = 100                  
 
    q = 7.25 * clay + 2.52             
    xwo = 0.33 * clay + 0.078          
    solidContent = bulkDensity / 2650
    porosity = 1 - solidContent
    gasPorosity = max(porosity - waterContent, 0.0)
    
    temperatureK = temperature + 273.16                     
    Lv = 45144 - 48 * temperature                 
    svp = 0.611 * exp(17.502 * temperature / (temperature + 240.97))
    slope = 17.502 * 240.97 * svp / (240.97 + temperature)**2.0
    Dv = 0.0000212 * (101.3 / atmPressure) * (temperatureK / 273.16)**1.75
    rhoair = 44.65 * (atmPressure / 101.3) * (273.16 / temperatureK)
    stcor = max(1 - svp / atmPressure, 0.3)

    thermalConductivitywater = 0.56 + 0.0018 * temperature
    if waterContent < 0.01 * xwo :
        wf = 0
    else:
        wf = 1 / (1 + (waterContent / xwo)**(-q))        
        
    thermalConductivitygas = (0.0242 + 0.00007 * temperature + 
                              wf * Lv * rhoair * Dv * slope / (atmPressure * stcor))
    gc = 1 - 2 * ga
    thermalConductivityfluid = (thermalConductivitygas + 
    (thermalConductivitywater - thermalConductivitygas) * (waterContent / porosity)**2.0)
    ka = (2 / (1 + (thermalConductivitygas / thermalConductivityfluid - 1) * ga) + 1 / 
         (1 + (thermalConductivitygas / thermalConductivityfluid - 1) * gc)) / 3
    kw = (2 / (1 + (thermalConductivitywater / thermalConductivityfluid - 1) * ga) + 
          1 / (1 + (thermalConductivitywater / thermalConductivityfluid - 1) * gc)) / 3
    ks = (2 / (1 + (thermalConductivitysolid / thermalConductivityfluid - 1) * ga) + 
          1 / (1 + (thermalConductivitysolid / thermalConductivityfluid - 1) * gc)) / 3
    
    thermalConductivity = ((kw * thermalConductivitywater * waterContent + ka * 
    thermalConductivitygas * gasPorosity + ks * thermalConductivitysolid * solidContent) 
                           / (kw * waterContent + ka * gasPorosity + ks * solidContent))
    return(thermalConductivity)       


def heatCapacity(bulkDensity, waterContent):
    return (2.4e6 * bulkDensity / 2650 + 4.18e6 * waterContent)       
