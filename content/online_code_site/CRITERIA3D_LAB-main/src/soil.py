# soil.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

import numpy as np
from math import fabs, log, sqrt
from dataStructures import *
import pandas as pd


class CSoilHorizon:
    upperDepth = NODATA     # [m]
    lowerDepth = NODATA     # [m]
    sand = NODATA           # [%]
    silt = NODATA           # [%]
    clay = NODATA           # [%]
    Campbell_he = NODATA    # [m] absolute value
    Campbell_b = NODATA     # [-]
    Campbell_n = NODATA     # [-]
    VG_he = NODATA          # [m] absolute value
    VG_alpha = NODATA       # [m-1]
    VG_n = NODATA           # [-]
    VG_m = NODATA           # [-]
    VG_Sc = NODATA          # [-]
    VG_thetaR = NODATA      # [m3 m-3]
    Mualem_L = NODATA       # [-]
    thetaS = NODATA         # [m3 m-3] water content at saturation
    Ks = NODATA             # [m s-1] soil water conductivity
    FC = NODATA             # [m3 m-3] water content at Field Capacity
    WP = NODATA             # [m3 m-3] water content at Wilting Point
    HYGR = NODATA           # [m3 m-3] water content at Hygroscopic moisture content


# global arrays
depth = np.array([], np.float64)
thickness = np.array([], np.float64)
horizons = []


def readHorizon(soilFileName):
    global horizons

    horizons.clear()
    soilDataFrame = pd.read_csv(soilFileName)

    i = 0
    for _, soilData in soilDataFrame.iterrows():
        horizon = CSoilHorizon()
        horizons.append(horizon)

        horizons[i].upperDepth = soilData["upper_depth"]
        # check upper depth
        if i > 0 and not horizons[i].upperDepth == horizons[i-1].lowerDepth:
            print("Wrong upper depth in layer: ", i+1)
            return False
        horizons[i].lowerDepth = soilData["lower_depth"]

        horizons[i].sand = soilData["sand"]
        horizons[i].silt = soilData["silt"]
        horizons[i].clay = soilData["clay"]

        horizons[i].Campbell_b = soilData["Campbell_b"]
        horizons[i].Campbell_n = 2.0 + (3.0 / horizons[i].Campbell_b)

        if soilData["Campbell_he"] < 0:
            horizons[i].Campbell_he = -soilData["Campbell_he"]
        else:
            horizons[i].Campbell_he = soilData["Campbell_he"]

        if soilData["VG_he"] < 0:
            horizons[i].VG_he = -soilData["VG_he"]
        else:
            horizons[i].VG_he = soilData["VG_he"]

        horizons[i].VG_alpha = soilData["VG_alpha"]
        horizons[i].VG_n = soilData["VG_n"]
        horizons[i].VG_m = 1. - (1. / horizons[i].VG_n)
        horizons[i].VG_Sc = pow(1. + pow(horizons[i].VG_alpha * fabs(horizons[i].VG_he),
                                         horizons[i].VG_n), -horizons[i].VG_m)

        horizons[i].VG_thetaR = soilData["thetaR"]
        horizons[i].thetaS = soilData["thetaS"]

        horizons[i].Ks = soilData["Ks"]
        horizons[i].Mualem_L = 0.5

        horizons[i].FC = getFieldCapacityWC(horizons[i])
        horizons[i].WP = getWiltingPointWC(horizons[i])
        horizons[i].HYGR = getHygroscopicWC(horizons[i])

        i += 1

    return True


def searchProgressionFactor(minThickness, maxThickness, maxThicknessDepth):
    if minThickness == maxThickness:
        return 1.0

    factor = 1.01
    bestError = 9999
    bestFactor = factor
    while factor <= 2.0:
        myThickness = minThickness
        currentDepth = minThickness * 0.5
        while myThickness < maxThickness:
            nextThickness = min(maxThickness, myThickness * factor)
            currentDepth += (myThickness + nextThickness) * 0.5
            myThickness = nextThickness
        error = fabs(currentDepth - maxThicknessDepth)
        if error < bestError:
            bestError = error
            bestFactor = factor
        factor += 0.01

    return bestFactor


# set depth and thickness of layers
def setLayers(totalDepth, minThickness, maxThickness, maxThicknessDepth):
    # search progression factor
    factor = searchProgressionFactor(minThickness, maxThickness, maxThicknessDepth)

    nrLayers = 1
    prevThickness = minThickness
    currentDepth = minThickness * 0.5
    while currentDepth < totalDepth:
        nextThickness = min(maxThickness, prevThickness * factor)
        currentDepth += (prevThickness + nextThickness) * 0.5
        prevThickness = nextThickness
        nrLayers += 1

    z = np.zeros(nrLayers, np.float64)
    thick = np.zeros(nrLayers, np.float64)
    z[0] = 0.0
    thick[0] = 0.0
    for i in range(1, nrLayers):
        top = z[i - 1] + thick[i - 1] * 0.5
        if i == 1:
            thick[i] = minThickness
        else:
            if i == (nrLayers - 1):
                thick[i] = totalDepth - top
            else:
                thick[i] = min(maxThickness, thick[i - 1] * factor)
        z[i] = top + thick[i] * 0.5
    return nrLayers, z, thick


def getHorizonIndex(myDepth):
    for i in range(len(horizons)):
        if horizons[i].upperDepth <= myDepth <= horizons[i].lowerDepth:
            return i
    return NODATA


def getVolumetricWaterContent(i):
    if C3DCells[i].isSurface:
        return NODATA
    curve = C3DParameters.waterRetentionCurve
    Se = C3DCells[i].Se
    index = C3DCells[i].horizonIndex
    return waterContent(curve, Se, horizons[index])


def getDegreeOfSaturation(i):
    if C3DCells[i].isSurface:
        if C3DCells[i].H > C3DCells[i].z:
            return 1.0
        else:
            return 0.0
    curve = C3DParameters.waterRetentionCurve
    index = C3DCells[i].horizonIndex
    signPsi = C3DCells[i].H - C3DCells[i].z
    return degreeOfSaturation(curve, signPsi, horizons[index])


def getHydraulicConductivity(i):
    if C3DCells[i].isSurface:
        return NODATA
    curve = C3DParameters.waterRetentionCurve
    index = C3DCells[i].horizonIndex
    return hydraulicConductivity(curve, C3DCells[i].Se, horizons[index])


# [m] air entry potential with sign
def airEntryPotential(curve, horizon):
    if curve == CAMPBELL:
        return -horizon.Campbell_he
    elif curve == IPPISCH_VG:
        return -horizon.VG_he
    else:
        return NODATA


def waterPotential(curve, Se, horizon):
    if curve == CAMPBELL:
        return -horizon.Campbell_he * Se ** (-horizon.Campbell_b)
    elif curve == IPPISCH_VG:
        return -(1. / horizon.VG_alpha) * ((1. / (Se * horizon.VG_Sc))
                                           ** (1. / horizon.VG_m) - 1.) ** (1. / horizon.VG_n)
    else:
        return NODATA


def waterContent(curve, Se, horizon):
    if curve == CAMPBELL:
        return Se * horizon.thetaS
    elif curve == IPPISCH_VG:
        return Se * (horizon.thetaS - horizon.VG_thetaR) + horizon.VG_thetaR
    else:
        return NODATA


def degreeOfSaturation(curve, signPsi, horizon):
    airEntry = airEntryPotential(curve, horizon)
    if signPsi >= airEntry:
        return 1.0

    Se = NODATA
    if curve == CAMPBELL:
        Se = pow(fabs(signPsi) / fabs(horizon.Campbell_he), -1. / horizon.Campbell_b)
    elif curve == IPPISCH_VG:
        Se = (1. / horizon.VG_Sc) * pow(1. + pow(horizon.VG_alpha
                                                 * fabs(signPsi), horizon.VG_n), -horizon.VG_m)
    return Se


def hydraulicConductivity(curve, Se, horizon):
    k = NODATA

    if curve == CAMPBELL:
        psi = horizon.Campbell_he * Se ** (-horizon.Campbell_b)
        k = horizon.Ks * (horizon.Campbell_he / psi) ** horizon.Campbell_n

    if curve == IPPISCH_VG:
        num = 1. - pow(1. - pow(Se * horizon.VG_Sc, 1. / horizon.VG_m), horizon.VG_m)
        den = 1. - pow(1. - pow(horizon.VG_Sc, 1. / horizon.VG_m), horizon.VG_m)
        k = horizon.Ks * pow(Se, horizon.Mualem_L) * pow((num / den), 2.)
    return k


def psiFromTheta(curve, theta, horizon):
    Se = SeFromTheta(curve, theta, horizon)
    return waterPotential(curve, Se, horizon)


def thetaFromPsi(curve, signPsi, horizon):
    Se = degreeOfSaturation(curve, signPsi, horizon)
    return waterContent(curve, Se, horizon)


def SeFromTheta(curve, theta, horizon):
    if theta >= horizon.thetaS:
        return 1.
    if curve == CAMPBELL:
        return theta / horizon.thetaS
    elif curve == IPPISCH_VG:
        return (theta - horizon.VG_thetaR) / (horizon.thetaS - horizon.VG_thetaR)
    else:
        return NODATA


def dTheta_dPsi(curve, signPsi, horizon):
    airEntry = airEntryPotential(curve, horizon)
    if signPsi > airEntry:
        return 0.0
    if curve == CAMPBELL:
        theta = horizon.thetaS * degreeOfSaturation(curve, signPsi, horizon)
        return -theta / (horizon.Campbell_b * signPsi)
    elif curve == IPPISCH_VG:
        dSe_dPsi = horizon.VG_alpha * horizon.VG_n * \
                   (horizon.VG_m * pow(1. + pow(horizon.VG_alpha * fabs(signPsi), horizon.VG_n), -(horizon.VG_m + 1.))
                    * pow(horizon.VG_alpha * fabs(signPsi), horizon.VG_n - 1.))
        dSe_dPsi *= (1. / horizon.VG_Sc)
        return dSe_dPsi * (horizon.thetaS - horizon.VG_thetaR)


def get_dTheta_dH(i):
    if C3DCells[i].isSurface:
        return NODATA
    curve = C3DParameters.waterRetentionCurve
    index = C3DCells[i].horizonIndex
    return dTheta_dH(curve, C3DCells[i].H0, C3DCells[i].H, C3DCells[i].z, horizons[index])


def dTheta_dH(curve, H0, H1, z, horizon):
    psi0 = H0 - z
    psi1 = H1 - z
    if fabs(psi1 - psi0) < EPSILON_METER:
        return dTheta_dPsi(curve, psi0, horizon)
    else:
        theta0 = thetaFromPsi(curve, psi0, horizon)
        theta1 = thetaFromPsi(curve, psi1, horizon)
        return (theta1 - theta0) / (psi1 - psi0)


def meanK(meanType, k1, k2):
    k = NODATA
    if meanType == LOGARITHMIC:
        if k1 != k2:
            k = (k1 - k2) / log(k1 / k2)
        else:
            k = k1
    elif meanType == HARMONIC:
        k = 2.0 / (1.0 / k1 + 1.0 / k2)
    elif meanType == GEOMETRIC:
        k = sqrt(k1 * k2)
    return k


# [m3 m-3] water content at field capacity
def getFieldCapacityWC(horizon):
    curve = C3DParameters.waterRetentionCurve
    fcMin = -10     # [kPa] clay < 20% : sandy soils
    fcMax = -33     # [kPa] clay > 50% : clay soils
    clayMin = 20    # [%]
    clayMax = 50    # [%]

    if horizon.clay < clayMin:
        fieldCapacity = fcMin
    elif horizon.clay >= clayMax:
        fieldCapacity = fcMax
    else:
        clayFactor = (horizon.clay - clayMin) / (clayMax - clayMin)
        fieldCapacity = fcMin + (fcMax - fcMin) * clayFactor

    FC = fieldCapacity / 9.81  # [m]
    return thetaFromPsi(curve, FC, horizon)


# [m3 m-3] water content at wilting point
def getWiltingPointWC(horizon):
    curve = C3DParameters.waterRetentionCurve
    WP = -1600      # [kPa]
    WP /= 9.81      # [m]
    return thetaFromPsi(curve, WP, horizon)


# [m3 m-3] water content at hygroscopic moisture
def getHygroscopicWC(horizon):
    curve = C3DParameters.waterRetentionCurve
    HH = -3100      # [kPa]
    HH /= 9.81      # [m]
    return thetaFromPsi(curve, HH, horizon)
