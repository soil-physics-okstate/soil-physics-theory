# criteria3D.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from math import fabs
import pandas as pd
import numpy as np
import time

from dataStructures import *
from PenmanMonteith import computeHourlyET0
import rectangularMesh
import waterBalance
import visual3D
import soil
import crop

if CYTHON:
    import solverCython as solver
else:
    import solver as solver

isRedraw = True


def setIsRedraw(value):
    global isRedraw
    isRedraw = value


def memoryAllocation(nrLayers, nrRectangles):
    C3DStructure.nrRectangles = nrRectangles
    C3DStructure.nrLayers = nrLayers
    C3DStructure.nrCells = nrLayers * nrRectangles
    solver.setCriteria3DArrays(C3DStructure.nrCells, C3DStructure.nrMaxLinks)
    C3DCells.clear()
    for i in range(C3DStructure.nrCells):
        C3DCells.append(Ccell())


def setCellGeometry(i, x, y, z, volume, area):
    C3DCells[i].x = x
    C3DCells[i].y = y
    C3DCells[i].z = z
    C3DCells[i].volume = volume
    C3DCells[i].area = area


def setCellProperties(i, isSurface, boundaryType):
    C3DCells[i].isSurface = isSurface
    C3DCells[i].boundary.type = boundaryType


def setBoundaryProperties(i, area, slope):
    C3DCells[i].boundary.area = area
    C3DCells[i].boundary.slope = slope


def getCellDistance(i, j):
    v1 = [C3DCells[i].x, C3DCells[i].y, C3DCells[i].z]
    v2 = [C3DCells[j].x, C3DCells[j].y, C3DCells[j].z]
    return rectangularMesh.distance3D(v1, v2)


# -----------------------------------------------------------
# direction:         UP, DOWN, LATERAL
# interfaceArea      [m^2]            
# -----------------------------------------------------------
def SetCellLink(i, linkIndex, direction, interfaceArea):
    if direction == UP:
        C3DCells[i].upLink.index = linkIndex
        C3DCells[i].upLink.area = interfaceArea
        C3DCells[i].upLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return OK
    elif direction == DOWN:
        C3DCells[i].downLink.index = linkIndex
        C3DCells[i].downLink.area = interfaceArea
        C3DCells[i].downLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return OK
    elif direction == LATERAL:
        for j in range(C3DStructure.nrLateralLinks):
            if C3DCells[i].lateralLink[j].index == NOLINK:
                C3DCells[i].lateralLink[j].index = linkIndex
                C3DCells[i].lateralLink[j].area = interfaceArea
                C3DCells[i].lateralLink[j].distance = getCellDistance(i, linkIndex)
                return OK
    else:
        return LINK_ERROR


def initializeMesh():
    print("Set cell properties...")
    waterTableDepth = NODATA

    for i in range(C3DStructure.nrRectangles):
        [x, y, z] = rectangularMesh.C3DRM[i].centroid
        dzy = C3DStructure.slopeY * y
        elevation = C3DStructure.elevation - dzy
        totalPotential = elevation + C3DParameters.initialWaterPotential
        if C3DParameters.isWaterTable:
            waterTableDepth = elevation + C3DParameters.waterTableDepth

        for layer in range(C3DStructure.nrLayers):
            index = i + C3DStructure.nrRectangles * layer
            elevation = z - soil.depth[layer]
            volume = float(rectangularMesh.C3DRM[i].area * soil.thickness[layer])
            setCellGeometry(index, x, y, elevation, volume, rectangularMesh.C3DRM[i].area)
            C3DCells[index].horizonIndex = soil.getHorizonIndex(soil.depth[layer])
            if layer == 0:
                # surface
                if rectangularMesh.C3DRM[i].isBoundary and C3DParameters.isSurfaceRunoff:
                    setCellProperties(index, True, BOUNDARY_RUNOFF)
                    setBoundaryProperties(index, rectangularMesh.C3DRM[i].boundarySide,
                                          rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    setCellProperties(index, True, BOUNDARY_NONE)

                setMatricPotential(index, max(C3DParameters.initialWaterPotential, 0.0))

            elif layer == (C3DStructure.nrLayers - 1):
                # last layer
                if C3DParameters.isWaterTable:
                    setCellProperties(index, False, BOUNDARY_PRESCRIBEDTOTALPOTENTIAL)
                elif C3DParameters.isFreeDrainage:
                    setCellProperties(index, False, BOUNDARY_FREEDRAINAGE)
                else:
                    setCellProperties(index, False, BOUNDARY_NONE)

                setTotalPotential(index, max(totalPotential, waterTableDepth))
            else:
                if rectangularMesh.C3DRM[i].isBoundary and C3DParameters.isFreeLateralDrainage:
                    setCellProperties(index, False, BOUNDARY_FREELATERALDRAINAGE)
                    setBoundaryProperties(index, rectangularMesh.C3DRM[i].boundarySide * soil.thickness[layer],
                                          rectangularMesh.C3DRM[i].boundarySlope)
                else:
                    setCellProperties(index, False, BOUNDARY_NONE)

                setTotalPotential(index, max(totalPotential, waterTableDepth))

    print("Set links...")
    for i in range(C3DStructure.nrRectangles):
        # UP
        for layer in range(1, C3DStructure.nrLayers):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i
            linkIndex = index - C3DStructure.nrRectangles
            SetCellLink(index, linkIndex, UP, exchangeArea)
        # LATERAL
        for neighbour in rectangularMesh.C3DRM[i].neighbours:
            if neighbour != NOLINK:
                linkSide = rectangularMesh.getAdjacentSide(i, neighbour)
                for layer in range(C3DStructure.nrLayers):
                    if layer == 0:
                        # surface: boundary length [m]
                        exchangeArea = linkSide
                    else:
                        # sub-surface: boundary area [m2]
                        exchangeArea = soil.thickness[layer] * linkSide
                    index = C3DStructure.nrRectangles * layer + i
                    linkIndex = C3DStructure.nrRectangles * layer + neighbour
                    SetCellLink(index, linkIndex, LATERAL, exchangeArea)
        # DOWN
        for layer in range(C3DStructure.nrLayers - 1):
            exchangeArea = rectangularMesh.C3DRM[i].area
            index = C3DStructure.nrRectangles * layer + i
            linkIndex = index + C3DStructure.nrRectangles
            SetCellLink(index, linkIndex, DOWN, exchangeArea)


def setTotalPotential(i, totalPotential):
    if C3DCells[i].isSurface:
        C3DCells[i].H = max(totalPotential, 0.0)
        C3DCells[i].Se = 1.
        C3DCells[i].k = soil.horizons[0].Ks
    else:
        C3DCells[i].H = totalPotential
        C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        C3DCells[i].k = soil.getHydraulicConductivity(i)
    C3DCells[i].H0 = C3DCells[i].H
    return OK


def setMatricPotential(i, signPsi):
    if C3DCells[i].isSurface:
        C3DCells[i].H = C3DCells[i].z + max(signPsi, 0.0)
        C3DCells[i].Se = 1.
        C3DCells[i].k = soil.horizons[0].Ks
    else:
        C3DCells[i].H = C3DCells[i].z + signPsi
        C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        C3DCells[i].k = soil.getHydraulicConductivity(i)
    C3DCells[i].H0 = C3DCells[i].H
    return OK


def getMatricPotential(i):
    if C3DCells[i].isSurface:
        return max(C3DCells[i].H - C3DCells[i].z, 0.0)
    else:
        return C3DCells[i].H - C3DCells[i].z


def initializeSinkSource(cellType):
    if cellType == ONLY_SURFACE:
        for i in range(C3DStructure.nrRectangles):
            C3DCells[i].sinkSource = 0
    else:
        for i in range(C3DStructure.nrCells):
            C3DCells[i].sinkSource = 0


# -----------------------------------------------------------
# set uniform rainfall rate
# rain            [mm]
# duration        [s]            
# -----------------------------------------------------------
def setRainfall(rain, duration):
    rate = (rain * 0.001) / duration  # [m s^-1]
    for i in range(C3DStructure.nrRectangles):
        area = C3DCells[i].area  # [m^2]
        C3DCells[i].sinkSource += rate * area  # [m^3 s^-1]


# -----------------------------------------------------------
# set drip irrigation
# irrigation      [l]
# duration        [s]            
# -----------------------------------------------------------
def setDripIrrigation(irrigation, duration):
    rate = irrigation / duration  # [l s^-1]

    for index in dripperIndices:
        C3DCells[index].sinkSource += rate * 0.001  # [m^3 s^-1]


# -----------------------------------------------------------
# compute water flow
# timeLength        [s]
# -----------------------------------------------------------
def computeWaterFlow(timeLength):
    currentTime = 0
    while currentTime < timeLength:
        residualTime = timeLength - currentTime
        deltaT = min(C3DParameters.currentDeltaT, residualTime)
        acceptedStep = False

        while not acceptedStep:
            if isRedraw:
                if visual3D.isPause and not visual3D.isComputeEquilibrium:
                    print("\nPress 'r' to run")
                    while visual3D.isPause:
                        time.sleep(0.00001)

            deltaT = min(C3DParameters.currentDeltaT, residualTime)
            # print("time step [s]: ", deltaT)
            # print("sink/source [l]:", format(waterBalance.sumSinkSource(deltaT) * 1000., ".5f"))

            acceptedStep = solver.computeStep(deltaT)
            if not acceptedStep:
                # restoreWater
                for i in range(C3DStructure.nrCells):
                    C3DCells[i].H = C3DCells[i].H0

        if isRedraw:
            visual3D.redraw()
        currentTime += deltaT


# -----------------------------------------------------------
# hourly computation model
# obsWeather    meteo data
# obsWater      precipitation [mm] and irrigation [l]
# transmissivity    normalized transmissivity [0-1]
# -----------------------------------------------------------
def computeOneHour(obsWeather, obsWater, transmissivity):
    currentDateTime = pd.to_datetime(obsWeather["timestamp"], unit='s')

    if not (np.isnan(obsWeather["air_temperature"])):
        airTemperature = obsWeather["air_temperature"]
    else:
        print("Missing data: airTemperature")
        return

    if not (np.isnan(obsWeather["solar_radiation"])):
        globalSWRadiation = obsWeather["solar_radiation"]
    else:
        print("Missing data: solar_radiation")
        return

    if not (np.isnan(obsWeather["air_humidity"])):
        airRelHumidity = obsWeather["air_humidity"]
    else:
        print("Missing data: air_humidity")
        return

    if not (np.isnan(obsWeather["wind_speed"])):
        windSpeed_10m = obsWeather["wind_speed"]
    else:
        print("Missing data: windspeed")
        windSpeed_10m = 2.0

    if C3DParameters.assignPrecipitation:
        if not (np.isnan(obsWater["precipitation"])):
            precipitation = obsWater["precipitation"]
        else:
            print("Missing data: precipitation")
            precipitation = 0
    else:
        precipitation = 0

    # evapotranspiration [mm m-2]
    ET0 = computeHourlyET0(C3DStructure.elevation, airTemperature, globalSWRadiation, airRelHumidity,
                           windSpeed_10m, transmissivity)

    initializeSinkSource(ALL)
    crop.setEvapotranspiration(currentDateTime, ET0)

    initializeSinkSource(ONLY_SURFACE)
    waterBalance.currentPrec = precipitation    # [mm hour-1]
    setRainfall(precipitation, 3600)

    if C3DParameters.assignIrrigation:
        if not (np.isnan(obsWater["irrigation"])):
            irrigation = obsWater["irrigation"]
        else:
            irrigation = 0
            print("Missing data: irrigation")

        waterBalance.currentIrr = irrigation    # [l hour-1]
        setDripIrrigation(irrigation, 3600)

    # reduce deltaT during water event
    if (waterBalance.currentIrr > 0) or (waterBalance.currentPrec > 0):
        if C3DParameters.currentDeltaT_max > 300:
            C3DParameters.currentDeltaT_max = 300
            C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT, C3DParameters.currentDeltaT_max)
    else:
        C3DParameters.currentDeltaT_max = C3DParameters.deltaT_max

    print(currentDateTime)
    computeWaterFlow(3600)


# -----------------------------------------------------------
# run until hydrostatic equilibrium
# -----------------------------------------------------------
def computeEquilibrium():
    initializeSinkSource(ALL)
    previousStorage = waterBalance.currentStep.waterStorage
    previousWaterFlow = waterBalance.allSimulation.waterFlow
    MBR = 999
    hour = 0
    print("hour:" + str(hour) + " water storage [m3]:" + format(previousStorage, ".6f"))
    while MBR > 1E-5:
        computeWaterFlow(3600)
        hour += 1
        currentStorage = waterBalance.currentStep.waterStorage
        deltaStorage = currentStorage - previousStorage

        if previousStorage == 0:
            MBR = abs(deltaStorage)
        else:
            MBR = abs(deltaStorage) / previousStorage

        hourFlow = waterBalance.allSimulation.waterFlow - previousWaterFlow
        print("hour:" + str(hour) + " water storage [m3]:" + format(currentStorage, ".6f")
               + " water flow [m3]:" + format(hourFlow, ".7f"))

        previousStorage = waterBalance.currentStep.waterStorage
        previousWaterFlow = waterBalance.allSimulation.waterFlow
