#PSP_boundary
from __future__ import division
import PSP_soil as soil
from PSP_public import *
import numpy as np

energy_balance = open('energy_balance.dat','w')
energy_balance.write('#Time \t LWUp \t SWDown \t LWDown \t Latent \t Sensible \t G \n')

class Cboundary:
    time = 0
    airTemperature = NODATA           
    precipitation = NODATA           
    relativeHumidity = NODATA
    windSpeed = NODATA
    globalRadiation = NODATA

def boundaryLayerConductance(windSpeed, Tair, Tsoil):
    z = 2.0                     
    cp = 29.3
    h = 0.01                   
    D = 0.77 * h                
    zm = 0.13 * h              
    zh = 0.2 * zm              
    psi_m = 0; psi_h = 0        
    vk = 0.4                   
    TsoilK = Tsoil + zeroKelvin     #[K]
    # molar density of the gas
    ro = 44.6 * (atmPressure / 101.3) * (293.15 / TsoilK)  
    # volumetric heat of air (= 1200 J/m^3*K at 20C e sea level)
    Ch = ro * cp                
    for i in range(3):
        #friction velocity        #[m/s]
        ustar = vk * windSpeed / (np.log((z - D + zm) / zm) + psi_m) 
        Kh = vk * ustar / (np.log((z - D + zh) / zh) + psi_h) 
        Sp = -vk * z * g * Kh * (Tair - Tsoil) / (Ch * TsoilK * np.power(ustar,3))
        if (Sp > 0):
            psi_h = 4.7 * Sp
            psi_m = psi_h
        else:
            psi_h = -2 * np.log((1 + np.sqrt(1 - 16 * Sp)) / 2)
            psi_m = 0.6 * psi_h
    return Kh

# soil resistance
def soilResistance(theta):
    return 10. * np.exp(0.3563 * (22. - (theta * 100.)))
 
# [kg/m2s] 
def evaporationFlux(psi, theta, Tsoil, Tair, rhAir, aerodynamicResistance):
    vapourConcAir = soil.vapourConcentration(Tair, rhAir/100.)  
    rhSoil = soil.relativeHumidity(psi, Tsoil + zeroKelvin)
    vapourConcSoil = soil.vapourConcentration(Tsoil, rhSoil) 
    dVapour = vapourConcSoil - vapourConcAir
    #print(vapourConcAir,rhSoil,vapourConcSoil, dVapour,aerodynamicResistance)
    return -(1.0 / (aerodynamicResistance + soilResistance(theta)) * dVapour) 

def dEvaporationFluxdPsi(psi, theta, Tsoil, Tair, aerodynamicResistance):
    TKelvinAir = Tair + zeroKelvin
    rhSoil = soil.relativeHumidity(psi, Tsoil + zeroKelvin)
    vapourConcSoil = soil.vapourConcentration(Tsoil, rhSoil)  
    return (1.0 / (aerodynamicResistance + soilResistance(theta)) * Mw/(R*TKelvinAir) * vapourConcSoil)
 
# [kg/m2s] 
def waterFlux(psi, theta, Tsoil, myBoundary, aerodynamicResistance):
    evapFlux = evaporationFlux(psi, theta, Tsoil, myBoundary.airTemperature, 
    myBoundary.relativeHumidity, aerodynamicResistance)
    precFlux = myBoundary.precipitation/3600
    return evapFlux + precFlux

def dWaterFluxdPsi(psi, theta, Tsoil, myBoundary, aerodynamicResistance):
    return (dEvaporationFluxdPsi(psi, theta, Tsoil, 
    myBoundary.airTemperature, aerodynamicResistance))

def thermalFlux(psi, theta, Tsoil, myBoundary, aerodynamicResistance, isWrite):
    TairK = myBoundary.airTemperature + zeroKelvin
    longWaveSoilTaylorAtmosphericPart =-(sigma*(TairK)**4 - 4.0*sigma*(TairK)**3 * myBoundary.airTemperature)
    sensibleHeatAtmosphericPart = 1200.0 / aerodynamicResistance * myBoundary.airTemperature  
    sensibleHeat = 1200.0 / aerodynamicResistance * (myBoundary.airTemperature-Tsoil)   
    evaporation = evaporationFlux(psi, theta, Tsoil, myBoundary.airTemperature, myBoundary.relativeHumidity, aerodynamicResistance)
    shortWaveAbsRadiation = (1 - albedo) * myBoundary.globalRadiation
    longWaveAbsRadiation = myBoundary.longWaveRadiation
    netAbsRadiation = shortWaveAbsRadiation + longWaveAbsRadiation 
    if isWrite:
        dTemperature = myBoundary.airTemperature - Tsoil
        sensibleHeat = (1200.0 / aerodynamicResistance * dTemperature) 
        longWaveSoil = -(sigma*(TairK)**4 -4.0*sigma*(TairK)**3 * dTemperature)
        G = longWaveAbsRadiation + shortWaveAbsRadiation + sensibleHeat + (evaporation*L) + longWaveSoil
        energy_balance.write(str(myBoundary.time) + '\t' 
                             + format(longWaveSoil, ".1f") 
                             + '\t' + format(shortWaveAbsRadiation, ".1f") 
                             + '\t' + format(longWaveAbsRadiation, ".1f") 
                             + '\t' + format(evaporation*L, ".1f") 
                             + '\t' + format(sensibleHeat, ".1f")
                             + '\t' + format(-G, ".1f") + '\n')
    return (longWaveSoilTaylorAtmosphericPart + netAbsRadiation + evaporation*L + sensibleHeatAtmosphericPart)