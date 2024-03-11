import numpy as np
from PSP_dataStructures import *
from PSP_readDataFile import readDataFile
import PSP_soil as soil
import PSP_balance as balance
import PSP_tin as tin
import PSP_criteria3D as Criteria3D
import PSP_visual3D as visual3D
from PSP_fileUtilities import loadState

 
def main():     
    print ("Load TIN...")
    vertexList, isFileOk = readDataFile("data/vertices.csv", 0, ",", False)
    if (not isFileOk): return 
    triangleList, isFileOk = readDataFile("data/triangles.csv", 0, ",", False) 
    if (not isFileOk): return 
    nrTriangles = len(triangleList)
    neighbourList, isFileOk = readDataFile("data/neighbours.csv", 0, ",", False) 
    if (not isFileOk): return 
    print ("Nr. of triangles:", nrTriangles)
    
    v = np.zeros((3, 3), float)
    for i in range(nrTriangles):
        for j in range(3):
            v[j] = vertexList[int(triangleList[i,j])]
        tin.C3DTIN.append(tin.Ctriangle(v))
        C3DStructure.totalArea += tin.C3DTIN[i].area
    print ("Total area [m^2]:", C3DStructure.totalArea)
    
    tin.header = tin.getHeader(tin.C3DTIN)
    
    print ("Set boundary...")
    for i in range(nrTriangles):
        tin.C3DTIN[i].isBoundary = False
        if (neighbourList[i,2] == NOLINK):
            tin.C3DTIN[i].isBoundary = True
            tin.boundaryProperties(tin.C3DTIN, i, neighbourList[i])
    
    # SOIL
    print ("Load soil...")
    soil.C3DSoil = soil.readHorizon("data/soil.txt", 1)
    if (C3DParameters.computeOnlySurface):
        totalDepth = 0
    else:
        totalDepth = soil.C3DSoil.lowerDepth
    print("Soil depth [m]:", totalDepth)
    
    nrLayers, soil.depth, soil.thickness = soil.setLayers(totalDepth, 
                     C3DParameters.minThickness, C3DParameters.maxThickness, 
                     C3DParameters.geometricFactor) 
    print("Nr. of layers:", nrLayers)
    
    # Initialize memory
    Criteria3D.memoryAllocation(nrLayers, nrTriangles)
    print("Nr. of cells: ", C3DStructure.nrCells)
    
    print("Set cell properties...")   
    for i in range(nrTriangles):
        for layer in range(nrLayers): 
            [x, y, z] = tin.C3DTIN[i].centroid 
            index = i + nrTriangles * layer
            elevation = z - soil.depth[layer]
            volume = float(tin.C3DTIN[i].area * soil.thickness[layer])
            Criteria3D.setCellGeometry(index, x, y, 
                                elevation, volume, tin.C3DTIN[i].area)
            if (layer == 0):
                # surface 
                if tin.C3DTIN[i].isBoundary:
                    Criteria3D.setCellProperties(index, True, BOUNDARY_RUNOFF)
                    Criteria3D.setBoundaryProperties(index, 
                                  tin.C3DTIN[i].boundarySide, tin.C3DTIN[i].boundarySlope)
                else:
                    Criteria3D.setCellProperties(index, True, BOUNDARY_NONE)
                Criteria3D.setMatricPotential(index, 0.0)
                
            elif (layer == (nrLayers-1)):
                # last layer
                if C3DParameters.isFreeDrainage:
                    Criteria3D.setCellProperties(index, False, BOUNDARY_FREEDRAINAGE)
                else:
                    Criteria3D.setCellProperties(index, False, BOUNDARY_NONE)
                    
                Criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)
                
            else:
                if tin.C3DTIN[i].isBoundary: 
                    Criteria3D.setCellProperties(index, False, BOUNDARY_FREELATERALDRAINAGE)
                    Criteria3D.setBoundaryProperties(index, tin.C3DTIN[i].boundarySide * soil.thickness[layer], 
                                                     tin.C3DTIN[i].boundarySlope)
                else:
                    Criteria3D.setCellProperties(index, False, BOUNDARY_NONE)
                    
                Criteria3D.setMatricPotential(index, C3DParameters.initialWaterPotential)
                 
    print("Set links...")   
    for i in range(nrTriangles): 
        # UP
        for layer in range(1, nrLayers):
            exchangeArea = tin.C3DTIN[i].area
            index = nrTriangles * layer + i 
            linkIndex = index - nrTriangles
            Criteria3D.SetCellLink(index, linkIndex, UP, exchangeArea)   
        # LATERAL
        for j in range(len(neighbourList[i])):
            neighbour = int(neighbourList[i,j])
            if (neighbour != NOLINK):
                linkSide = tin.getAdjacentSide(i, neighbour, vertexList, triangleList)
                for layer in range(nrLayers): 
                    if (layer == 0):
                        #surface: boundary length [m]
                        exchangeArea = linkSide
                    else:
                        #sub-surface: boundary area [m2]
                        exchangeArea = soil.thickness[layer] * linkSide
                    index = nrTriangles * layer + i 
                    linkIndex = nrTriangles * layer + neighbour
                    Criteria3D.SetCellLink(index, linkIndex, LATERAL, exchangeArea)
        # DOWN
        for layer in range(nrLayers-1):
            exchangeArea = tin.C3DTIN[i].area
            index = nrTriangles * layer + i 
            linkIndex = index + nrTriangles
            Criteria3D.SetCellLink(index, linkIndex, DOWN, exchangeArea)
            
    # LOAD INITIAL STATE - comment if you dont't have one
    if (not C3DParameters.computeOnlySurface): 
        print ("Load initial state...")
        loadState("data/state_0.csv")
    
    balance.initializeBalance()
    print("Initial water storage [m^3]:", format(balance.currentStep.waterStorage, ".3f"))
        
    print("Read precipitation data...")
    precFileName = "data/precipitation.txt"
    # TIME LENGHT 
    # change it if your observed data are different (for example: hourly)
    timeLength = 15 * 60         # [s]
    C3DParameters.deltaT_max = timeLength
    print("Time lenght [s]:", timeLength)
    data, isFileOk = readDataFile(precFileName, 1, "\t", False)
    if (not isFileOk):
        print("Error! Wrong precipitation file.") 
        return
    prec = data[:,1]
    nrObsPrec = len(prec)
    print("Total simulation time [s]:", nrObsPrec * timeLength)
    
    visual3D.initialize(1200)
    visual3D.isPause = True
    
    # main cycle
    for i in range(nrObsPrec):
        balance.currentPrec = prec[i] / timeLength * 3600   #[mm/hour]
        Criteria3D.setRainfall(prec[i], timeLength)
        Criteria3D.compute(timeLength)
    
    print ("\nEnd simulation.")      
main()
