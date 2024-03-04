#PSP_longWaveRadiation
from __future__ import division
import numpy as np
from PSP_public import *
from math import pi

def vaporConcentrationAir(airT, relativeHumidity):  
    c_vsat = 0.611*1.e3*np.exp(17.27*airT/(airT+237.3)) * Mw/(R*(airT+273.15))
    c_v = relativeHumidity * c_vsat
    return c_v

def atmEmissivity(doy, dailyRad, airT, relativeHumidity):
    latRad = latitude * (pi / 180.0)
    sin_sD = (0.3985 * np.sin(4.869 + 0.0172*doy + 0.03345 * np.sin(6.224 + 0.0172*doy)))  
    cos_sD = np.sqrt(1.0 - sin_sD * sin_sD) 
    h_s = np.arccos(-np.tan(latRad) * (sin_sD/cos_sD))
    potentialRadiation = (117.5e6 * (h_s*np.sin(latRad) * sin_sD
			+np.cos(latRad) * cos_sD * np.sin(h_s)) / 3.1415)
    T_t = dailyRad / potentialRadiation
    c_1 = 2.33 - 3.33 * T_t
    if(c_1 < 0): 
	    c_1 = 0.
    elif(c_1 > 1):
	    c_1 = 1.
    c_va = vaporConcentrationAir(airT, relativeHumidity) * 1.e3  
    epsilon_a = 0.58 * np.power(c_va, 1./7.)
    emissivity = (1. - 0.84*c_1) * epsilon_a + 0.84 * c_1
    return emissivity

def longWaveRadiationFromWeather(firstDoy, airT, relativeHumidity, globalRadiation):
    nrHours = len(airT)
    nrDays = int(nrHours / 24.)
    longWaveRadiation = np.zeros(nrHours)
    for d in range(nrDays):
            dailyRad = sum(globalRadiation[d*24:(d+1)*24]) * 3600.0
            for j in range(24):
                i = d*24+j
                doy = firstDoy + d
                longWaveRadiation[i] = (atmEmissivity(doy, dailyRad, airT[i], relativeHumidity[i]/100.)
                                    *sigma *np.power(airT[i]+273.15,4.))
    return longWaveRadiation
