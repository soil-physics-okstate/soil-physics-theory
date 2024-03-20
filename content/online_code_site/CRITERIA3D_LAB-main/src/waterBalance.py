# waterBalance.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from math import fabs
from dataStructures import *
from soil import getVolumetricWaterContent


class C3DBalance:
    waterStorage = NODATA       # [m3]
    waterFlow = NODATA          # [m3]
    MBE = NODATA                # [m3] Mass Balance Error
    MBR = NODATA                # [-] Mass Balance Ratio


totalTime = 0.0
currentPrec = 0.0
currentIrr = 0.0
MBRMultiply: float = 1.0
maxCourant = 0.0
bestMBR = NODATA
nrMBRWrong = 0
forceExit = False
currentStep = C3DBalance()
previousStep = C3DBalance()
allSimulation = C3DBalance()


def doubleTimeStep():
    global MBRMultiply
    if (C3DParameters.currentDeltaT == C3DParameters.deltaT_min) and (MBRMultiply > 1.0):
        decMBRThreshold()
    else:
        C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT * 2.0,
                                          C3DParameters.currentDeltaT_max)


def halveTimeStep():
    if C3DParameters.currentDeltaT == C3DParameters.deltaT_min:
        if C3DParameters.MBRThreshold < 1E-2:
            incMBRThreshold()
        else:
            return False
    else:
        C3DParameters.currentDeltaT = max(C3DParameters.currentDeltaT * 0.5,
                                          C3DParameters.deltaT_min)
    return True


def incMBRThreshold():
    global MBRMultiply
    MBRMultiply *= 2.0
    C3DParameters.MBRThreshold *= 2.0


def decMBRThreshold():
    global MBRMultiply
    if MBRMultiply > 1.0:
        MBRMultiply *= 0.5
        C3DParameters.MBRThreshold *= 0.5


def initializeBalance():
    global totalTime

    totalTime = 0.0
    storage = getWaterStorage()
    currentStep.waterStorage = storage
    previousStep.waterStorage = storage
    allSimulation.waterStorage = storage
    previousStep.waterFlow = 0.0
    currentStep.waterFlow = 0.0
    allSimulation.waterFlow = 0.0
    currentStep.MBR = 0.0
    currentStep.MBE = 0.0
    allSimulation.MBE = 0


def updateStorage():
    storage = getWaterStorage()
    currentStep.waterStorage = storage
    previousStep.waterStorage = storage
    allSimulation.waterStorage = storage


def updateBalance(deltaT):
    global totalTime
    totalTime += deltaT
    previousStep.waterStorage = currentStep.waterStorage
    previousStep.waterFlow = currentStep.waterFlow
    allSimulation.waterFlow += currentStep.waterFlow
    allSimulation.MBE += currentStep.MBE


def getWaterStorage():
    waterStorage = 0.0
    for i in range(C3DStructure.nrCells):
        if C3DCells[i].isSurface:
            if abs(C3DCells[i].H - C3DCells[i].z) > EPSILON:
                waterStorage += (C3DCells[i].H - C3DCells[i].z) * C3DCells[i].area
        else:
            waterStorage += (getVolumetricWaterContent(i) * C3DCells[i].volume)
    return waterStorage


def sumBoundaryFlow(deltaT):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if C3DCells[i].boundary.type != BOUNDARY_NONE:
            if C3DCells[i].boundary.flow != NODATA:
                mySum += C3DCells[i].boundary.flow * deltaT
    return mySum


def sumSinkSource(deltaT):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if C3DCells[i].sinkSource != NODATA:
            mySum += C3DCells[i].sinkSource * deltaT
    return mySum


# [m3]
def sumWaterFlow(deltaT, isAbsoluteValue):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if C3DCells[i].flow != NODATA:
            if isAbsoluteValue:
                mySum += abs(C3DCells[i].flow * deltaT)
            else:
                mySum += C3DCells[i].flow * deltaT
    return mySum


def computeBalanceError(deltaT):
    currentStep.waterStorage = getWaterStorage()
    currentStep.waterFlow = sumWaterFlow(deltaT, False)
    currentStep.MBE = currentStep.waterStorage - (previousStep.waterStorage + currentStep.waterFlow)
    if previousStep.waterStorage > 0:
        currentStep.MBR = fabs(currentStep.MBE) / previousStep.waterStorage
    else:
        currentStep.MBR = fabs(currentStep.MBE)

    # print ("Mass Balance Error [l]:", format(currentStep.MBE * 1000,".5f"))
    # print("Mass Balance Ratio:", format(currentStep.MBR, ".5f"))


def waterBalance(deltaT, approximation):
    global forceExit, bestMBR, nrMBRWrong
    computeBalanceError(deltaT)

    if approximation == 1:
        bestMBR = currentStep.MBR
        nrMBRWrong = 0
        forceExit = False

    # case 1: step accepted
    if currentStep.MBR < C3DParameters.MBRThreshold:
        updateBalance(deltaT)
        if approximation < 3 and maxCourant < 0.3 and currentStep.MBR < (C3DParameters.MBRThreshold * 0.5) \
                and C3DParameters.currentDeltaT < C3DParameters.currentDeltaT_max:
            # print("Good MBR!")
            doubleTimeStep()
        return True

    # case 2: continue with next approximation
    if approximation == 1 or currentStep.MBR < bestMBR:
        bestMBR = currentStep.MBR
    else:
        nrMBRWrong += 1

    # case 3: decrease time step (or increase threshold)
    isLastApprox = (approximation == C3DParameters.maxApproximationsNr)
    if isLastApprox or nrMBRWrong > 0:
        if halveTimeStep():
            forceExit = True
        else:
            # accept error
            updateBalance(deltaT)
            return True

    return False
