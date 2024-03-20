#PSP_soil.py
from math import fabs, log, sqrt
import numpy as np
from PSP_readDataFile import readDataFile
from PSP_dataStructures import *


class CsoilHorizon:
    upperDepth = NODATA        # [m]
    lowerDepth = NODATA        # [m]
    Campbell_he = NODATA       # [J kg^-1] 
    Campbell_b = NODATA        # [-] 
    Campbell_n = NODATA        # [-] 
    VG_he = NODATA             # [J kg^-1] 
    VGmod_alpha = NODATA       # [kg J^-1]
    VGmod_n = NODATA           # [-]
    VGmod_m = NODATA           # [-]
    VG_Sc = NODATA             # [-] 
    VG_thetaR = NODATA         # [m^3 m^-3] 
    Mualem_L = NODATA          # [-] 
    thetaS = NODATA            # [m^3 m^-3] 
    Ks = NODATA                # [kg s m^-3]

#global
depth = np.array([], np.float64)
thickness = np.array([], np.float64)
C3DSoil = CsoilHorizon()

def readHorizon(soilFileName, i):
    A, isFileOk = readDataFile(soilFileName,1,',', False)
    if ((not isFileOk) or (len(A[0]) < 10)):
        print("warning: wrong soil file.")
        return False
    horizon = CsoilHorizon()
    i -= 1
    horizon.upperDepth = A[i,0]
    horizon.lowerDepth = A[i,1]
    horizon.Campbell_he = -A[i,2]
    horizon.Campbell_b = A[i,3]
    horizon.Campbell_n = 2.0 + (3.0 / horizon.Campbell_b)
    horizon.VG_he = -A[i,4]
    horizon.VG_alpha = A[i,5]
    horizon.VG_n = A[i,6]
    horizon.VG_m =  1. - (1. / horizon.VG_n)
    horizon.VG_Sc = pow(1. + pow(horizon.VG_alpha * fabs(horizon.VG_he), 
                                horizon.VG_n), -horizon.VG_m)
    horizon.VG_thetaR = A[i,7]
    horizon.thetaS = A[i,8]
    horizon.Ks = A[i,9]  
    horizon.Mualem_L = 0.5
    return horizon

# set depth and thickness of layers
def setLayers(totalDepth, minThickness, maxThickness, factor):
    nrLayers = 1
    prevThickness = minThickness
    depth = minThickness * 0.5
    while (depth < totalDepth):
        nextThickness = min(maxThickness, prevThickness * factor)
        depth = depth + (prevThickness + nextThickness) * 0.5
        prevThickness = nextThickness
        nrLayers += 1
        
    z = np.zeros(nrLayers, np.float64)
    thickness = np.zeros(nrLayers, np.float64)    
    z[0] = 0.0
    thickness[0] = 0.0
    for i in range(1, nrLayers):
        if (i == 1):
            thickness[i] = minThickness
        else:
            thickness[i] = min(maxThickness, thickness[i-1] * factor)
        z[i] = z[i-1] + (thickness[i-1] + thickness[i]) * 0.5
    return(nrLayers, z, thickness)

def getVolumetricWaterContent(i):
    if C3DCells[i].isSurface: return NODATA
    curve = C3DParameters.waterRetentionCurve
    Se = C3DCells[i].Se
    return waterContent(curve, Se)

def getDegreeOfSaturation(i):
    if C3DCells[i].isSurface:
        if  (C3DCells[i].H > C3DCells[i].z): return 1.0
        else: return 0.0
    curve = C3DParameters.waterRetentionCurve
    signPsi = C3DCells[i].H - C3DCells[i].z
    return degreeOfSaturation(curve, signPsi)

def getHydraulicConductivity(i):
    if C3DCells[i].isSurface: return NODATA
    curve = C3DParameters.waterRetentionCurve
    return hydraulicConductivity(curve, C3DCells[i].Se)

def airEntryPotential(curve): 
    if (curve == CAMPBELL):
        return C3DSoil.Campbell_he
    elif (curve == IPPISCH_VG):
        return C3DSoil.VG_he
    else: return NODATA
    
def waterPotential(curve, Se):
    if (curve == CAMPBELL):
        return C3DSoil.Campbell_he * Se**(-C3DSoil.Campbell_b)
    elif (curve == IPPISCH_VG):
        return -(1./C3DSoil.VG_alpha)*((1./(Se*C3DSoil.VG_Sc))
                 **(1./C3DSoil.VG_m)-1.)**(1./C3DSoil.VG_n)
    else: return NODATA
    
def waterContent(curve, Se):
    if (curve == CAMPBELL):
        return Se * C3DSoil.thetaS 
    elif (curve == IPPISCH_VG):
        return (Se * (C3DSoil.thetaS - C3DSoil.VG_thetaR) + C3DSoil.VG_thetaR)
    else: return NODATA
    
def degreeOfSaturation(curve, signPsi):
    airEntry = airEntryPotential(curve)
    if (signPsi >= airEntry): return 1.0
    
    Se = NODATA
    if (curve == CAMPBELL):
        Se = pow(signPsi / C3DSoil.Campbell_he, -1. / C3DSoil.Campbell_b)
    elif (curve == IPPISCH_VG):
        Se = (1./C3DSoil.VG_Sc) * pow(1. + pow(C3DSoil.VG_alpha 
                  * fabs(signPsi), C3DSoil.VG_n), -C3DSoil.VG_m)
    return Se

def hydraulicConductivity(curve, Se): 
    k = NODATA    
    if (curve == CAMPBELL):
        psi = C3DSoil.Campbell_he * Se**(-C3DSoil.Campbell_b)
        k = C3DSoil.Ks * (C3DSoil.Campbell_he / psi)**C3DSoil.Campbell_n  
    if (curve == IPPISCH_VG):
        num   = 1. - pow(1. - pow(Se * C3DSoil.VG_Sc, 1./ C3DSoil.VG_m), C3DSoil.VG_m)
        denom = 1. - pow(1. - pow(C3DSoil.VG_Sc, 1./ C3DSoil.VG_m), C3DSoil.VG_m)
        k = C3DSoil.Ks * pow(Se, C3DSoil.Mualem_L) * pow((num / denom), 2.)
    return(k)
    
def psiFromTheta(curve, theta):
    Se = SeFromTheta(curve, theta)
    return waterPotential(curve, Se)

def thetaFromPsi(curve, signPsi):
    Se = degreeOfSaturation(curve, signPsi)
    return waterContent(curve, Se)
    
def SeFromTheta(curve, theta):
    if (theta >= C3DSoil.thetaS): return(1.)
    if (curve == CAMPBELL):
        return theta / C3DSoil.thetaS
    elif (curve == IPPISCH_VG):
        return (theta - C3DSoil.VG_thetaR) / (C3DSoil.thetaS - C3DSoil.VG_thetaR)
    else: return NODATA 

def dTheta_dPsi(curve, signPsi):
    airEntry = airEntryPotential(curve)
    if (signPsi > airEntry): return 0.0
    if (curve == CAMPBELL):
        theta = C3DSoil.thetaS * degreeOfSaturation(curve, signPsi) 
        return -theta / (C3DSoil.Campbell_b * signPsi)
    elif (curve == IPPISCH_VG):
        dSe_dpsi = C3DSoil.VG_alpha * C3DSoil.VG_n * (C3DSoil.VG_m 
                * pow(1. + pow(C3DSoil.VG_alpha * fabs(signPsi), C3DSoil.VG_n), 
                -(C3DSoil.VG_m + 1.)) * pow(C3DSoil.VG_alpha * fabs(signPsi), 
                                                        C3DSoil.VG_n - 1.))      
        dSe_dpsi *= (1. / C3DSoil.VG_Sc)
        return dSe_dpsi * (C3DSoil.thetaS - C3DSoil.VG_thetaR)
    
def getdTheta_dH(i):
    if C3DCells[i].isSurface: return NODATA
    curve = C3DParameters.waterRetentionCurve
    return dTheta_dH(curve, C3DCells[i].H0, C3DCells[i].H, C3DCells[i].z)

def dTheta_dH(curve, H0, H1, z): 
    psi0 = H0 - z
    psi1 = H1 - z
    if (fabs(psi1-psi0) < EPSILON_METER):
        return dTheta_dPsi(curve, psi0)
    else:
        theta0 = thetaFromPsi(curve, psi0)
        theta1 = thetaFromPsi(curve, psi1)
        return ((theta1 - theta0) / (psi1 - psi0))
    
def meanK(meanType, k1, k2):
    if (meanType == LOGARITHMIC):
        if (k1 != k2):
            k = (k1-k2) / log(k1/k2)
        else:
            k = k1
    elif (meanType == HARMONIC): 
        k = 2.0 / (1.0 / k1 + 1.0 / k2)
    elif (meanType == GEOMETRIC): 
        k = sqrt(k1 * k2)
    return k

