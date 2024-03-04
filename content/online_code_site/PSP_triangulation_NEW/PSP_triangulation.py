#PSP_triangulation.py
from __future__ import print_function, division
from PSP_gis import *
import PSP_refinement as refinement
from PSP_RamerDouglasPeucker import getDtmBoundary
from PSP_Delaunay import firstTriangulation 

import sys
if sys.version_info < (3, 0):
    import PSP_visual3D_27 as visual3D          #2.7
else:
    import PSP_visual3D as visual3D             #3.x                


def searchPosition(x, sortedList, first, last): 
    if x <= sortedList[first][0]:
        return first
    elif x > sortedList[last][0]:
        return (last+1)
    elif (last - first) < 2:
        return last
    else:
        m = int((first+last)/2)
        if (x <= sortedList[m][0]):
            return(searchPosition(x, sortedList, first, m))
        else:
            return(searchPosition(x, sortedList, m, last))


def sortPointListX(pointList):
    sortedList = [pointList[0]]
    
    for i in range(1, len(pointList)):
        x = pointList[i][0]
        index = searchPosition(x, sortedList, 0, len(sortedList) - 1)
        sortedList.insert(index, pointList[i])
        
    return(sortedList)


def initializeFlagMatrix(header, dtm, pointList):
    flagMatrix = np.zeros((header.nrCols, header.nrRows), bool)
    
    for col in range(header.nrCols):
        for row in range(header.nrRows):
            if dtm[col][row] == header.flag:
                flagMatrix[col, row] = False
            else:
                flagMatrix[col, row] = True
                  
    for i in range(len(pointList)): 
        p = pointList[i]
        col, row = getColRowFromXY(header, p[0], p[1])   
        if not isOutOfGridColRow(header, col, row):
            flagMatrix[col, row] = False
            
    return flagMatrix   


def addPartitionPoints(pointList, header, dtm, intervalStep, flagMatrix):
    print ("Add partition points...")
    col = -1
    isLastCol = False
    while not isLastCol:
        row = -1
        isLastRow = False
        while not isLastRow:
            c = col + intervalStep
            r = row + intervalStep
            if isTrue(header, flagMatrix, c, r): 
                point = getPointFromColRow(header, dtm, c, r)
                pointList.append(point)
                flagMatrix[c, r] = False
            row += intervalStep
            if row >= (header.nrRows-1): isLastRow = True
        col += intervalStep
        if col >= (header.nrCols-1): isLastCol = True
        

def triangulate(header, dtm, compressionRate, angleThreshold, initialPartition):             
    # [points] average distance between two points of the initial partition
    partitionStep = round(sqrt(sqrt(header.nrPoints)))
    
    # [m] 
    boundaryStep = header.cellSize * sqrt(compressionRate)
    
    # [m] maximum height difference for refinement   
    thresholdZ = (header.dz / 100.0) * sqrt(compressionRate)
    
    # [m] maximum perpendicular distance for boundary   
    thresholdRamer = thresholdZ * sqrt(compressionRate)
     
    # [m^2] minimum area of triangles
    minArea = header.cellSize**2 * compressionRate
    
    pointList = getDtmBoundary(header, dtm, thresholdRamer, boundaryStep)
    flagMatrix = initializeFlagMatrix(header, dtm, pointList)
     
    if initialPartition:    
        addPartitionPoints(pointList, header, dtm, partitionStep, flagMatrix)
        
    pointList = sortPointListX(pointList)
    pointList, triangleList = firstTriangulation(pointList, header, dtm)
    
    visual3D.drawAllTriangles(triangleList, header)
    
    pointList, triangleList, flagMatrix = refinement.refinementZ(pointList, triangleList, 
                        header, dtm, thresholdZ, minArea, flagMatrix)
    
    visual3D.drawAllTriangles(triangleList, header)
    
    pointList, triangleList, flagMatrix = refinement.refinementAngle(pointList, triangleList, 
                        header, dtm, angleThreshold, minArea, flagMatrix)
    
    pointList = sortPointListX(pointList)  
    return (pointList, triangleList)

