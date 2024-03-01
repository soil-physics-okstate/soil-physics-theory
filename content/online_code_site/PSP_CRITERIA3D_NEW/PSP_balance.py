#PSP_balance.py
from math import fabs
from PSP_dataStructures import *
from PSP_public import EPSILON_METER, C3DParameters
from PSP_soil import getVolumetricWaterContent

class C3DBalance:
    waterStorage = NODATA
    waterFlow = NODATA
    MBE = NODATA
    MBR = NODATA

totalTime = 0.0
currentPrec = 0.0
MBRMultiply = 1.0
maxCourant = 0.0
bestMBR = NODATA
forceExit = False
currentStep = C3DBalance()
previousStep = C3DBalance()
allSimulation = C3DBalance()
    
def doubleTimeStep():
    C3DParameters.currentDeltaT = min(C3DParameters.currentDeltaT * 2.0, 
                                      C3DParameters.deltaT_max)
    
def halveTimeStep():
    C3DParameters.currentDeltaT = max(C3DParameters.currentDeltaT * 0.5, 
                                      C3DParameters.deltaT_min)
    
def incMBRThreshold():
    global MBRMultiply
    MBRMultiply *= 2.0
    C3DParameters.MBRThreshold *= 2.0
    
def decMBRThreshold():
    global MBRMultiply
    if (MBRMultiply > 1.0):
        MBRMultiply *= 0.5
        C3DParameters.MBRThreshold *= 0.5
        
def initializeBalance():
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
        if (C3DCells[i].isSurface):
            if (C3DCells[i].H > C3DCells[i].z):
                waterStorage += (C3DCells[i].H - C3DCells[i].z) * C3DCells[i].area
        else:
            waterStorage += (getVolumetricWaterContent(i) * C3DCells[i].volume)
    return waterStorage
                    
def sumBoundaryFlow(deltaT):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if (C3DCells[i].boundary.type != BOUNDARY_NONE):
            if (C3DCells[i].boundary.flow != NODATA):
                mySum += C3DCells[i].boundary.flow * deltaT
    return (mySum)

def sumSinkSource(deltaT):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if (C3DCells[i].sinkSource != NODATA):
            mySum += C3DCells[i].sinkSource * deltaT
    return (mySum)

def sumWaterFlow(deltaT, isAbsoluteValue):
    mySum = 0.0
    for i in range(C3DStructure.nrCells):
        if (C3DCells[i].flow != NODATA):
            if isAbsoluteValue:
                mySum += fabs(C3DCells[i].flow * deltaT)
            else:
                mySum += C3DCells[i].flow * deltaT
    return (mySum)
                                    
def computeBalanceError(deltaT):
    currentStep.waterStorage = getWaterStorage()
    currentStep.waterFlow = sumWaterFlow(deltaT, False)
    deltaStorage = currentStep.waterStorage - previousStep.waterStorage
    currentStep.MBE = deltaStorage - currentStep.waterFlow
    
    sumFlow = sumWaterFlow(deltaT, True)
    minimumFlow = C3DStructure.totalArea * EPSILON_METER * (deltaT / 3600.0)
    if (sumFlow < minimumFlow):
        currentStep.MBR = fabs(currentStep.MBE) / minimumFlow
    else:
        currentStep.MBR = fabs(currentStep.MBE) / sumFlow
    print ("Mass Balance Error [m^3]:", format(currentStep.MBE,".3f"))
    print ("Mass Balance Ratio:", format(currentStep.MBR,".5f"))
    
def waterBalance(deltaT, approximation):
    global bestMBR, forceExit
    if (approximation == 1): bestMBR = 100.0
    computeBalanceError(deltaT)
    isLastApprox = (approximation == C3DParameters.maxApproximationsNr)
    
    forceExit = False
    # case 1: error < tolerance
    if currentStep.MBR <= C3DParameters.MBRThreshold:
        updateBalance(deltaT)
        if ((approximation < 3) and (maxCourant < 0.5) 
        and (currentStep.MBR < (C3DParameters.MBRThreshold * 0.5))):
                print("Good MBR!")
                if (deltaT > (C3DParameters.deltaT_min * 2)):
                    decMBRThreshold()
                doubleTimeStep()
        return True
    #case 2: error decreases
    if (currentStep.MBR < bestMBR):
        bestMBR = currentStep.MBR
        if isLastApprox:
            if (deltaT == C3DParameters.deltaT_min):
                updateBalance(deltaT)
                return True
            else:
                halveTimeStep()
                forceExit = True
                return False 
        return False
    # case 3: error increases
    if (deltaT > C3DParameters.deltaT_min):
        print("Solution is not convergent: decrease time step")
        halveTimeStep()
        forceExit = True
        return False
    else:
        print("Solution is not convergent: increase error tolerance")
        incMBRThreshold()
        return False
        
        