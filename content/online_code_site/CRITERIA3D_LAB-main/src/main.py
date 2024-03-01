# main.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

import pandas as pd
import os
import time

from dataStructures import *
import soil
import waterBalance
import criteria3D
import visual3D
import exportUtils
import importUtils
import crop
from transmissivity import computeNormTransmissivity


def main():
    # path
    print(os.getcwd())
    projectPath = os.path.join("..\\data", "test2D")
    settingsFolder = os.path.join(projectPath, "settings")
    weatherFolder = os.path.join(projectPath, "meteo")
    obsDataFolder = os.path.join(projectPath, "obs_data")
    stateFolder = os.path.join(projectPath, "state")
    outputFolder = os.path.join(projectPath, "output")

    print("Read model settings...")
    modelSettings = os.path.join(settingsFolder, "settings.ini")
    if not importUtils.readModelParameters(modelSettings):
        return

    print("Read field settings...")
    fieldSettings = os.path.join(settingsFolder, "field.ini")
    if not importUtils.readFieldParameters(fieldSettings):
        return

    print("read soil properties...")
    soilSettings = os.path.join(settingsFolder, "soil.csv")
    if not soil.readHorizon(soilSettings):
        return
    if C3DStructure.gridDepth > soil.horizons[len(soil.horizons)-1].lowerDepth:
        print("Wrong soil properties: lower depth is < field.depth")
        return
    C3DStructure.nrLayers, soil.depth, soil.thickness = soil.setLayers(C3DStructure.gridDepth,
                                                                       C3DParameters.minThickness,
                                                                       C3DParameters.maxThickness,
                                                                       C3DParameters.maxThicknessAt)
    print("Nr. of layers:", C3DStructure.nrLayers)

    criteria3D.memoryAllocation(C3DStructure.nrLayers, C3DStructure.nrRectangles)
    print("Nr. of cells: ", C3DStructure.nrCells)

    # initialize crop
    if C3DParameters.computeTranspiration:
        cropSettingsFilename = os.path.join(settingsFolder, "crop.ini")
        if os.path.exists(cropSettingsFilename):
            print("Read crop settings...")
            if not importUtils.readCropParameters(cropSettingsFilename):
                return
            crop.initializeCrop()
        else:
            print("WARNING: crop settings file does not exist!")
            print("*** The transpiration process will be deactivated.")
            C3DParameters.computeTranspiration = False
            crop.setNoCrop()
    else:
        crop.setNoCrop()

    print("Initialize mesh...")
    criteria3D.initializeMesh()

    waterBalance.initializeBalance()
    print("Initial water storage [m3]:", format(waterBalance.currentStep.waterStorage, ".5f"))

    # read unified meteo input file
    print("Read weather and irrigation data...")
    wholeData = importUtils.readMeteoData(weatherFolder)
    wholeData["time"] = pd.to_datetime(wholeData["timestamp"], infer_datetime_format=True)

    print("Total simulation time [hours]:", len(wholeData))

    # split the dataframe in two
    weatherData = wholeData[["timestamp", "air_humidity", "solar_radiation", "air_temperature", "wind_speed"]]
    waterData = wholeData[["timestamp", "precipitation", "irrigation"]]

    # initialize export
    exportUtils.createExportFile(outputFolder)

    obsTmpFileName = os.path.join(stateFolder, "obsWP.csv")
    if C3DParameters.isPeriodicAssimilation or C3DParameters.isFirstAssimilation:
        obsWaterPotentialFileName = os.path.join(obsDataFolder, "waterPotential.csv")
        if os.path.exists(obsWaterPotentialFileName):
            print("Read observed water potential...")
            obsWaterPotential = pd.read_csv(obsWaterPotentialFileName)
        else:
            print("WARNING: observed water potential file does not exist!")
            print("*** The assimilation procedure will be deactivated.")
            C3DParameters.isPeriodicAssimilation = False
            C3DParameters.isFirstAssimilation = False

    # first assimilation
    weatherIndex = 0
    if C3DParameters.isFirstAssimilation:
        print("Assimilate observed water potential (first hour)...")
        obsWeather = weatherData.loc[weatherIndex]
        timestamp = obsWeather["timestamp"]
        if not importUtils.extractObsWaterPotential(obsWaterPotential, timestamp, obsTmpFileName):
            return
        importUtils.assimilateObsWaterPotential(obsTmpFileName)
        criteria3D.setIsRedraw(False)
        criteria3D.computeWaterFlow(3600)
        importUtils.assimilateObsWaterPotential(obsTmpFileName)

    criteria3D.setIsRedraw(C3DParameters.isVisual)
    if C3DParameters.isVisual:
        visual3D.initialize(1200)
        visual3D.isPause = True
        # wait for start (press 'r')
        while visual3D.isPause:
            time.sleep(0.00001)
            # check for equilibrium
            if visual3D.isComputeEquilibrium:
                criteria3D.computeEquilibrium()
                visual3D.isComputeEquilibrium = False

    print("Start...")
    waterBalance.initializeBalance()
    if C3DParameters.isVisual:
        visual3D.redraw()

    # main cycle
    currentIndex = 1
    lastIndex = min(len(weatherData), len(waterData))
    while weatherIndex < lastIndex:
        # compute
        obsWeather = weatherData.loc[weatherIndex]
        obsWater = waterData.loc[weatherIndex]
        transmissivity = computeNormTransmissivity(weatherData, weatherIndex, C3DStructure.latitude, C3DStructure.longitude)
        criteria3D.computeOneHour(obsWeather, obsWater, transmissivity)

        # assimilation
        if C3DParameters.isPeriodicAssimilation and (currentIndex % C3DParameters.assimilationInterval) == 0:
            print("Assimilate observed water potential...")
            if not importUtils.extractObsWaterPotential(obsWaterPotential, obsWeather["timestamp"], obsTmpFileName):
                return
            importUtils.assimilateObsWaterPotential(obsTmpFileName)

        # save output
        exportUtils.takeScreenshot(obsWeather["timestamp"])

        weatherIndex += 1
        currentIndex += 1

    visual3D.isPause = True
    print("\nEnd simulation.\n")


main()
