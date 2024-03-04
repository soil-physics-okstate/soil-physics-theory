#PSP_criteria3D.py
from math import fabs
from PSP_dataStructures import *
import PSP_soil as soil
from PSP_balance import sumSinkSource
from PSP_tin import distance3D
import PSP_visual3D as visual3D
import time

CYTHON = True
if CYTHON:
    import PSP_solverCython as solver
else:
    import PSP_solver as solver


def memoryAllocation(nrLayers, nrTriangles):
    C3DStructure.nrTriangles = nrTriangles
    C3DStructure.nrLayers = nrLayers
    nrCells = nrLayers * nrTriangles
    C3DStructure.nrCells = nrCells
    solver.setCriteria3DArrays(nrCells, C3DStructure.nrMaxLinks)
    for i in range(nrCells): 
        C3DCells.append(Ccell())

def setCellGeometry(i, x, y, z, volume, area):
    C3DCells[i].x = x;
    C3DCells[i].y = y;
    C3DCells[i].z = z;
    C3DCells[i].volume = volume;
    C3DCells[i].area = area;

def setCellProperties(i, isSurface, boundaryType):
    C3DCells[i].isSurface = isSurface
    C3DCells[i].boundary.type = boundaryType

def setBoundaryProperties(i, area, slope):
    C3DCells[i].boundary.area = area
    C3DCells[i].boundary.slope = slope

def getCellDistance(i, j):
    v1 = [C3DCells[i].x, C3DCells[i].y, C3DCells[i].z]
    v2 = [C3DCells[j].x, C3DCells[j].y, C3DCells[j].z]
    return distance3D(v1, v2)

#-----------------------------------------------------------
# direction:         UP, DOWN, LATERAL
# interfaceArea      [m^2]            
#-----------------------------------------------------------
def SetCellLink(i, linkIndex, direction, interfaceArea):
    if (direction == UP):
        C3DCells[i].upLink.index = linkIndex
        C3DCells[i].upLink.area = interfaceArea
        C3DCells[i].upLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return(OK)
    elif(direction == DOWN):
        C3DCells[i].downLink.index = linkIndex
        C3DCells[i].downLink.area = interfaceArea
        C3DCells[i].downLink.distance = fabs(C3DCells[i].z - C3DCells[linkIndex].z)
        return(OK)
    elif (direction == LATERAL):
        for j in range(C3DStructure.nrLateralLinks):
            if (C3DCells[i].lateralLink[j].index == NOLINK):     
                C3DCells[i].lateralLink[j].index = linkIndex
                C3DCells[i].lateralLink[j].area = interfaceArea 
                C3DCells[i].lateralLink[j].distance = getCellDistance(i, linkIndex)
                return(OK)
    else:
        return(LINK_ERROR)

def setMatricPotential (i, signPsi):
    if (C3DCells[i].isSurface):
        C3DCells[i].H = C3DCells[i].z + max(signPsi, 0.0)
        C3DCells[i].Se = 1.
        C3DCells[i].k = soil.C3DSoil.Ks
    else: 
        C3DCells[i].H = C3DCells[i].z + signPsi
        C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        C3DCells[i].k = soil.getHydraulicConductivity(i)
    C3DCells[i].H0 = C3DCells[i].H
    return(OK)
                           
#-----------------------------------------------------------
# set uniform rainfall rate
# rain            [mm]
# duration        [s]            
#-----------------------------------------------------------
def setRainfall(rain, duration):   
    rate = (rain * 0.001) / duration                    #[m s^-1]
    for i in range(C3DStructure.nrTriangles):
        area = C3DCells[i].area                         #[m^2]
        C3DCells[i].sinkSource = rate * area            #[m^3 s^-1]
        
def restoreWater():
    for i in range(C3DStructure.nrCells):
        C3DCells[i].H = C3DCells[i].H0
        
# timeLength        [s]          
def compute(timeLength):  
    currentTime = 0
    while (currentTime < timeLength):  
        residualTime = timeLength - currentTime
        acceptedStep = False
        while (not acceptedStep):
            if visual3D.isPause:
                print ("\nPress 'r' to run")
            while visual3D.isPause:
                time.sleep(0.00001)
                
            deltaT = min(C3DParameters.currentDeltaT, residualTime)
            print ("\ntime step [s]: ", deltaT)
            print ("MBR threshold [-]: ", C3DParameters.MBRThreshold)
            print ("sink/source [m^3]:", format(sumSinkSource(deltaT),".3f")) 
             
            acceptedStep = solver.computeStep(deltaT)          
            if not acceptedStep: 
                restoreWater()
                
        visual3D.redraw()  
        currentTime += deltaT
