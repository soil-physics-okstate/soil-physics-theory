# PenmanMonteith.py
# ----------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ----------------------------------------------------------
#
# Penman Monteith equation for potential evapotranspiration
# References:
# Allen et al., 1998, Crop Evapotranspiration – 
# Guidelines for Computing Crop Water Requirements

import math

# [W m-2 K-4] Stefan-Boltzmann constant
STEFAN_BOLTZMANN = 5.670373E-8
# [Pa] standard atmospheric pressure at sea level
P0 = 101300.
# [K m-1] constant lapse rate of moist air
LAPSE_RATE_MOIST_AIR = 0.0065
# [J kg-1 K-1] specific gas constant for dry air
R_DRY_AIR = 287.058
# [J kg-1 K-1] specific heat at constant pressure
CP = 1013.
# [m s-2] gravity acceleration
GRAVITY = 9.80665
# [K]
ZERO_CELSIUS = 273.15
# [K] temperature at reference pressure level (P0)
TEMP_P0 = 20. + ZERO_CELSIUS
# [-] albedo of reference crop (grass)
ALBEDO_CROP_REFERENCE = 0.23
# [-] ratio molecular weight of water vapour/dry air
RATIO_WATER_VD = 0.622


# pressure [Pa]
# height [m]
def pressureFromAltitude(height):
    return P0 * math.pow(1. + height * LAPSE_RATE_MOIST_AIR / TEMP_P0,
                         - GRAVITY / (LAPSE_RATE_MOIST_AIR * R_DRY_AIR))


# saturation vapor pressure [Pa]
# airTemperature [degC] 
def saturationVaporPressure(airTemperature):
    return 611 * math.exp(17.502 * airTemperature / (airTemperature + 240.97))


# return reference evapotranspiration (mm)
# height               elevation above mean sea level (meters)
# airTemperature       air temperature (C)
# globalSWRadiation    global Short Wave radiation (W m-2)
# airRelHumidity       air relative humidity (%)
# windSpeed_10m        wind speed at 10 meters (m s-1)
# normTransmissivity   normalized transmissivity [0-1]
def computeHourlyET0(height, airTemperature, globalSWRadiation, airRelHumidity,
                     windSpeed_10m, normTransmissivity):
    # air temperature [Kelvin]  
    airTempKelvin = airTemperature + ZERO_CELSIUS
    # wind speed at 2 meters [m s-1]
    windSpeed_2m = windSpeed_10m * 0.748
    # barometric pressure [kPa]
    pressure = pressureFromAltitude(height) / 1000.
    # saturation vapor pressure [kPa] 
    satVapPressure = saturationVaporPressure(airTemperature) / 1000.
    # current vapor pressure [kPa] 
    vaporPressure = satVapPressure * (airRelHumidity / 100.)
    # net emissivity of the surface [-]
    emissivity = 0.34 - 0.14 * math.sqrt(vaporPressure)
    # cloudiness factor for long wave radiation [-]
    cloudFactor = max(0, 1.35 * min(normTransmissivity, 1) - 0.35)
    # net long wave radiation [J m-2 s-1]
    netLWRadiation = cloudFactor * emissivity * STEFAN_BOLTZMANN * math.pow(airTempKelvin, 4.)
    # net radiation [J m-2 s-1] 
    netRadiation = (1. - ALBEDO_CROP_REFERENCE) * globalSWRadiation - netLWRadiation

    # from [W m-2] to [J m-2 h-1]
    netRadiation = netRadiation * 3600.

    # values for grass
    # g   soil heat flux density [J m-2 h-1]  
    # Cd  bulk surface resistance and aerodynamic resistance coefficient  
    if netRadiation > 0:
        g = 0.1 * netRadiation
        Cd = 0.24
    else:
        g = 0.5 * netRadiation
        Cd = 0.96

    # slope of saturation vapor pressure curve [kPa K]
    slope = 4098. * satVapPressure / (airTempKelvin * airTempKelvin)
    # latent heat of vaporization [J kg-1]
    latentHeatVap = 2501000. - 2369.2 * airTemperature
    # psychrometric instrument constant [kPa K-1] 
    psychro = CP * pressure / (RATIO_WATER_VD * latentHeatVap)

    denominator = slope + psychro * (1. + Cd * windSpeed_2m)
    firstTerm = slope * (netRadiation - g) / (latentHeatVap * denominator)
    secondTerm = (psychro * (37. / airTempKelvin) * windSpeed_2m * (satVapPressure - vaporPressure)) / denominator

    return max(firstTerm + secondTerm, 0.)
