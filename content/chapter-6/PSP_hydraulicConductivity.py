#PSP_hydraulicConductivity.py
import numpy as np
from math import exp

def hydraulicConductivity(silt,clay,bulkDensity):
    psi = np.logspace(0, 3, 100) # logarithmic sequence of psi values
    
    def geoMeanDiam(silt,clay):
        sand = 1-silt-clay # sand mass fraction (g/g)
        gmd = exp(-1.96*clay+2.3*silt+5.76*sand) # geometric mean diameter (micrometers); Eq. 2.36
        return gmd
        
    def campbellParameters(gmd):
        psi_e = 0.61*np.log(gmd)-3.9 # air entry potential (J/kg); Eq. 5.36
        b = 8.25-1.26*np.log(gmd) # Campbell exponent (-); Eq.5.37
        return psi_e, b
    
    def computeK(bulkDensity, psi_e, b):
        porosity = 1-bulkDensity/2650; # porosity
        theta_s = porosity # saturated water content
        K_s = 0.07*(theta_s*(1-(-psi_e/33)**(1/b)))**4 #Eq. 6.34, cm/s
        K = np.piecewise(-psi, [-psi < psi_e, -psi >= psi_e], [lambda x: K_s*(psi_e/x)**(2+3/b), K_s]) # hydraulic conductivity (cm/s)
        return K
    
    gmd = geoMeanDiam(silt,clay)
    psi_e, b = campbellParameters(gmd)
    K = computeK(bulkDensity, psi_e, b)
    return K, psi



