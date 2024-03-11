#PSP_soil
from __future__ import division
from PSP_readDataFile import readDataFile
from PSP_public import *
import numpy as np

CAMPBELL = 1
RESTRICTED_VG = 2
IPPISCH_VG = 3
VAN_GENUCHTEN = 4

CELL_CENT_FIN_VOL = 1
NEWTON_RAPHSON_MP = 2
NEWTON_RAPHSON_MFP = 3

LOGARITHMIC = 0
HARMONIC = 1
GEOMETRIC = 2

class Csoil:
    upperDepth = NODATA        
    lowerDepth = NODATA        
    Campbell_he = NODATA       
    Campbell_b = NODATA        
    CampbellMFP_he = NODATA    
    VG_alpha = NODATA          
    VG_n = NODATA              
    VG_m = NODATA
    VG_he = NODATA             
    VG_alpha_mod = NODATA      
    VG_n_mod = NODATA          
    VG_m_mod = NODATA
    VG_Sc = NODATA             
    VG_thetaR = NODATA         
    Mualem_L = NODATA          
    thetaS = NODATA            
    Ks = NODATA                
    
def readSoil(soilFileName):
    soil = Csoil()
    A, isFileOk = readDataFile(soilFileName, 1, ',', False)
    if ((not isFileOk) or (len(A[0]) < 12)):
        return False, soil
    
    soil.upperDepth = A[0,0]
    soil.lowerDepth = A[0,1]
    soil.Campbell_he = A[0,2]
    soil.Campbell_b = A[0,3]
    soil.Campbell_n = 2.0 + (3.0 / soil.Campbell_b)
    soil.VG_he = A[0,4]
    soil.VG_alpha = A[0,5]
    soil.VG_n = A[0,6]
    soil.VG_m =  1. - (1. / soil.VG_n)
    soil.VG_alpha_mod = A[0,7]
    soil.VG_n_mod = A[0,8]
    soil.VG_m_mod =  1. - (1. / soil.VG_n_mod)
    soil.VG_Sc = ((1. + 
    (soil.VG_alpha_mod * 
    abs(soil.VG_he))**soil.VG_n_mod)**(-soil.VG_m_mod))
    soil.VG_thetaR = A[0,9]
    soil.thetaS = A[0,10]
    soil.Ks = A[0,11]
    soil.Mualem_L = 0.5
    soil.CampbellMFP_he = (soil.Ks * 
    soil.Campbell_he / (1.0 - soil.Campbell_n)) 
    return True, soil

def airEntryPotential(funcType, soil): 
    if (funcType == CAMPBELL):
        return(soil.Campbell_he)
    elif (funcType == IPPISCH_VG):
        return(soil.VG_he)
    elif (funcType == RESTRICTED_VG):
        return(0)
    else:
        return(NODATA)
     
def waterPotential(funcType, soil, theta):
    psi = NODATA
    Se = SeFromTheta(funcType, soil, theta)
    if (funcType == RESTRICTED_VG):
        psi = (-(1./soil.VG_alpha)*
        ((1./Se)**(1./soil.VG_m) - 1.)**(1./soil.VG_n))
    elif (funcType == IPPISCH_VG):
        psi = (-((1./soil.VG_alpha_mod)*
        ((1./(Se*soil.VG_Sc))
         **(1./soil.VG_m_mod)-1.)**(1./soil.VG_n_mod)))
    elif (funcType == CAMPBELL):
        psi = soil.Campbell_he * Se**(-soil.Campbell_b)
    return(psi)
    
def SeFromTheta(funcType, soil, theta):
    check = np.greater_equal(theta,soil.thetaS)
    if (funcType == CAMPBELL):
        Se = theta / soil.thetaS
    else:
        Se = ((theta - soil.VG_thetaR) 
        / (soil.thetaS - soil.VG_thetaR))
    return (1.*check + Se*np.logical_not(check))

def thetaFromSe(funcType, soil, Se):
    if (funcType == RESTRICTED_VG) or (funcType == IPPISCH_VG):
        theta = (Se * (soil.thetaS 
        - soil.VG_thetaR) + soil.VG_thetaR)
    elif (funcType == CAMPBELL):
        return(Se * soil.thetaS) 
    return(theta)

def degreeOfSaturation(funcType, soil, psi):
    if (funcType == IPPISCH_VG):
        Se = ((1./soil.VG_Sc) * pow(1.+pow(soil.VG_alpha_mod 
        * np.maximum(psi,0.), soil.VG_n_mod), -soil.VG_m_mod))
        return Se 
    elif (funcType == RESTRICTED_VG):
        Se = (1 / pow(1 + 
        pow(soil.VG_alpha * np.maximum(-psi,0.), 
        soil.VG_n), soil.VG_m))  
        return Se 
    elif (funcType == CAMPBELL):
        Se = (pow(np.minimum(psi,soil.Campbell_he) 
        / soil.Campbell_he, -1. / soil.Campbell_b))
        return Se

def thetaFromPsi(funcType, soil, psi):
    Se = degreeOfSaturation(funcType, soil, psi)
    theta = thetaFromSe(funcType, soil, Se)
    return(theta)
          
def dTheta_dPsi(funcType, soil, psi):
    airEntry = airEntryPotential(funcType, soil)
    if (funcType == RESTRICTED_VG):
        check = np.greater(psi,airEntry)
        dSe_dpsi = soil.VG_alpha * soil.VG_n * (soil.VG_m 
        * pow(1. + pow(soil.VG_alpha * abs(psi), soil.VG_n), 
        -(soil.VG_m + 1.)) * pow(soil.VG_alpha * abs(psi), soil.VG_n - 1.))      
        return (0*check + dSe_dpsi * (soil.thetaS 
        - soil.VG_thetaR)*np.logical_not(check))
    elif (funcType == IPPISCH_VG):
        check = np.greater(psi,airEntry)
        dSe_dpsi = soil.VG_alpha_mod * soil.VG_n_mod * (soil.VG_m_mod 
        * pow(1. + pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod), 
        -(soil.VG_m_mod + 1.)) * 
        pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod - 1.))      
        dSe_dpsi *= (1. / soil.VG_Sc)
        return (0*check + dSe_dpsi * (soil.thetaS 
        - soil.VG_thetaR)*np.logical_not(check))
    elif (funcType == CAMPBELL):
        check = np.greater(psi,airEntry)
        theta = soil.thetaS * degreeOfSaturation(funcType, soil, psi) 
        return (0*check + -theta 
        / (soil.Campbell_b * psi)*np.logical_not(check))

def hydraulicConductivityFromPsi(funcType, soil, psi):
    if (funcType == RESTRICTED_VG):
        psi = abs(psi)
        num = (1. - pow(soil.VG_alpha * psi, soil.VG_m * soil.VG_n)
        *pow(1. + pow(soil.VG_alpha*psi, soil.VG_n), -soil.VG_m))**2
        denom = (pow(1. + pow(soil.VG_alpha*psi, soil.VG_n), 
        soil.VG_m * soil.Mualem_L))
        k = soil.Ks * (num / denom)
    elif (funcType == IPPISCH_VG):
        k = NODATA
    elif (funcType == CAMPBELL):
        k = soil.Ks * (soil.Campbell_he / psi)**soil.Campbell_n
    return(k)

def relativeHumidity(psi, TKelvin):
    return np.exp(psi * Mw/(R*TKelvin))

def vapourDiffusivity(mySoil, theta, TKelvin, atmPressure):
    binaryDiffCoeff = (waterVapourDiff0 * (101.3 / atmPressure) 
    * (TKelvin / 273.15)**1.75)
    gasPorosity = (mySoil.thetaS - theta)
    bg = 0.9            
    mg = 2.3
    return (binaryDiffCoeff * bg * gasPorosity**mg)

def satVapourPressure(TCelsius): 
    return 0.611 * np.exp(17.502 * TCelsius / (TCelsius + 240.97))

def satVapourConcentration(TCelsius): 
    svp = satVapourPressure(TCelsius)
    TKelvin = TCelsius + zeroKelvin
    return svp * 1.e3 * Mw/(R*TKelvin)

def vapourPressure(TCelsius, rh):
    return satVapourPressure(TCelsius) * rh                 

def vapourConcentration(TCelsius, rh):
    return satVapourConcentration(TCelsius) * rh      

def slopeSatVapourConcentration(TCelsius, atmPressure,satVapConc):
    return 4098.0 * satVapConc / (TCelsius+237.3)**2              
    
def vapourConductivity(mySoil, psi, theta, TCelsius, atmPressure):
    TKelvin = TCelsius + zeroKelvin
    rh = relativeHumidity(psi, TKelvin)
    vapDiff = vapourDiffusivity(mySoil, theta, TKelvin, atmPressure)
    satVapConc = satVapourConcentration(TCelsius)
    return vapDiff * satVapConc * rh * Mw/(R*TKelvin)

def computeVapour(mySoil, psi, theta, rh, satVapConc):
    gasPorosity = (mySoil.thetaS - theta)
    return vapourConcentration(rh, satVapConc) * gasPorosity

def heatCapacity(mySoil, waterContent):
    solidContent = 1.0 - mySoil.thetaS
    return (2.4e6 * solidContent + 4.18e6 * waterContent) 

def thermalConductivity(mySoil, waterContent, TCelsius):
    ga = 0.088                         
    thermalConductivitysolid = 2.5      
 
    q = 7.25 * clay + 2.52              
    xwo = 0.33 * clay + 0.078          
    
    porosity = mySoil.thetaS
    solidContent = 1 - porosity
    gasPorosity = np.maximum(porosity - waterContent, 0.0)
    
    TKelvin = TCelsius + zeroKelvin           
    Lv = 45144 - 48 * TCelsius                
    svp = satVapourPressure(TCelsius)
    slope = 17.502 * 240.97 * svp / (240.97 + TCelsius)**2.0
    Dv = 0.0000212 * (101.3 / atmPressure) * (TKelvin / 273.16)**1.75
    rhoair = 44.65 * (atmPressure / 101.3) * (273.16 / TKelvin)
    stcor = np.maximum(1 - svp / atmPressure, 0.3)

    thermalConductivityWater = 0.56 + 0.0018 * TCelsius
    #empirical weighting function D[0,1]
    check = np.less(waterContent,0.01*xwo)
    wf = 0*check + (1 / (1 + (waterContent / xwo)**(-q)))*np.logical_not(check)   
        
    thermalConductivityGas = (0.0242 + 0.00007 * TCelsius + 
    wf * Lv * rhoair * Dv * slope / (atmPressure * stcor))
    gc = 1 - 2 * ga
    thermalConductivityFluid = (thermalConductivityGas + 
    (thermalConductivityWater - thermalConductivityGas) 
    * (waterContent / porosity)**2.0)
    ka = (2 / (1 + (thermalConductivityGas / thermalConductivityFluid - 1) * ga) + 1 / 
    (1 + (thermalConductivityGas / thermalConductivityFluid - 1) * gc)) / 3
    kw = (2 / (1 + (thermalConductivityWater / thermalConductivityFluid - 1) * ga) + 
    1 / (1 + (thermalConductivityWater / thermalConductivityFluid - 1) * gc)) / 3
    ks = (2 / (1 + (thermalConductivitysolid / thermalConductivityFluid - 1) * ga) + 
    1 / (1 + (thermalConductivitysolid / thermalConductivityFluid - 1) * gc)) / 3
    
    thermalConductivity = ((kw * thermalConductivityWater * waterContent + ka * 
    thermalConductivityGas * gasPorosity + ks * thermalConductivitysolid * solidContent) 
                           / (kw * waterContent + ka * gasPorosity + ks * solidContent))
    return(thermalConductivity)  

class PropertiesForMassBalance:
    def __init__(self,TCelsius,psi,mySoil):
        k = hydraulicConductivityFromPsi(CAMPBELL, mySoil, psi)
        theta = thetaFromPsi(CAMPBELL, mySoil, psi)
        TKelvin = TCelsius + zeroKelvin
        rh = relativeHumidity(psi, TKelvin)
        vapourDiff = vapourDiffusivity(mySoil, theta, TKelvin, atmPressure)
        satVapConc = satVapourConcentration(TCelsius)
        slope = slopeSatVapourConcentration(TCelsius, atmPressure,satVapConc)
        #vapour = computeVapour(mySoil, psi, theta, rh, satVapConc)
        thermalCond = thermalConductivity(mySoil,theta, TCelsius)
        
        self.K11 = np.append(0,k + vapourDiff*satVapConc*rh*Mw/(R*TKelvin))
        self.gravityWaterFlow = np.append(0,-k*g)                                  
        self.K21 = np.append(0,vapourDiff*rh*slope)                               
        self.K12 = np.append(0,L*vapourDiff*satVapConc*rh*Mw/(R*TKelvin))         
        self.K22 = np.append(0,thermalCond + L*vapourDiff*rh*slope)
       
