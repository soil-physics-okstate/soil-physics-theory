#PSP_TTwaterContent.py
from __future__ import division
from math import sqrt, exp

SPEED_OF_LIGHT = 299792458      # [m s^-1]                       
airPermittivity = 1.00058986    # [-]
waterDensity = 997              # [kg m^-3]

def getLiquidPermittivity(temperature):
    deltaT = temperature - 25.
    return(78.54 * (1-4.579E-03 * deltaT))
    
def getBulkPermittivity(probleLenght, travelTime, Vp):
    return(((SPEED_OF_LIGHT * Vp * travelTime) / (2. * probleLenght))**2)

def getWaterContentTopp(bulkPermittivity):
    return(-5.3E-02 + 2.92E-02 * bulkPermittivity - 5.5E-04 * bulkPermittivity**2
           + 4.3E-06 * bulkPermittivity**3)
    
def getWaterContentMalicki(bulkPermittivity, bulkDensity):
    bulkDensity /= 1000.
    return((sqrt(bulkPermittivity) - 0.819 - 0.168*bulkDensity - 0.159*bulkDensity**2)
            / (7.17 + 1.18*bulkDensity))
    
def getWaterContentMixModel(bulkPermittivity, bulkDensity, 
                            solidPermittivity, liquidPermittivity, alpha):
    porosity = 1. - bulkDensity/2650.
    numerator = bulkPermittivity**alpha - ((1. - porosity) * solidPermittivity**alpha 
                                           + porosity * airPermittivity**alpha)
    denominator =  liquidPermittivity**alpha - airPermittivity**alpha 
    return(numerator/denominator)

def getDensityJung(bulkPermittivity, v1, vf, c1, d1, f1):
    numerator = v1/vf
    denominator = c1 + d1*(bulkPermittivity-1) - c1*exp(-f1*(bulkPermittivity-1))
    return(waterDensity * numerator/denominator)

def getDensityCurioni(bulkPermittivity, v1, vf, a, b, c):
    vr = v1/vf
    numerator = waterDensity*vr
    denominator = a+b*(v1*sqrt(bulkPermittivity))**c
    return (numerator/denominator)
    
