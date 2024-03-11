#PSP_Delaunay.py
import numpy as np
from copy import copy
import PSP_triangle as triangle
from PSP_gis import NODATA


def insertVertexClockwise(newVertex, center, vertexList, angleList):
    #is first point
    v = copy(newVertex)
    if (len(vertexList) == 0):
        vertexList.append(v)
        angleList.append(-2)
        return()
    #check duplicate
    for i in range(len(vertexList)):
        if ((v[0] == vertexList[i][0])
        and (v[1] == vertexList[i][1])):
            return()
    #compute false angle (no trigonometric functions) [-2,2]    
    angle = triangle.getPseudoAngle(v, center, vertexList[0])
    #check position
    i = 1
    while ((i < len(angleList)) and (angle > angleList[i])):
        i += 1
    vertexList.insert(i, v)
    angleList.insert(i, angle)
    
    
def Delaunay(triangleList, newPoint):
    deleteList = []
    vertexList = []
    angleList = []
    newTriangles = []
    # erase a triangle if the point is inside his circumcircle.
    for i in range(len(triangleList)): 
        myTriangle = triangleList[i]    
        if triangle.isInCircumCircle(newPoint, myTriangle):
            deleteList.append(i)
            for j in range(3):
                insertVertexClockwise(myTriangle.v[j], 
                            newPoint, vertexList, angleList)                          
    # create new triangles
    nrVertices = len(vertexList) 
    for i in range(nrVertices):
        v = np.zeros((3, 3), float)
        v[0] = newPoint
        v[1] = vertexList[i]
        v[2] = vertexList[(i+1) % nrVertices] 
        myTriangle = triangle.Ctriangle(v)
        if not myTriangle.circle.isCorrect: return False
        newTriangles.append(myTriangle)
        
    for i in range(len(deleteList)-1, -1, -1):
        triangleList.pop(deleteList[i])
    for i in range(len(newTriangles)):
        triangleList.append(newTriangles[i])                   
    return True


def firstTriangulation(pointList, header, dtm): 
    print("First triangulation...") 
    step = header.cellSize * (header.nrCols + header.nrRows) / 50
    xMin = header.xllCorner - step
    xMax = header.xllCorner + header.nrCols * header.cellSize + step
    yMin = header.yllCorner - step
    yMax = header.yllCorner + header.nrRows * header.cellSize + step
    
    A = np.array([xMin, yMin, NODATA])
    B = np.array([xMin, yMax, NODATA])
    C = np.array([xMax, yMax, NODATA])
    D = np.array([xMax, yMin, NODATA])
    
    # boundary: domain divided in two triangles
    triangleList = []
    triangleList.append(triangle.Ctriangle(np.array([A,B,C])))
    triangleList.append(triangle.Ctriangle(np.array([C,D,A])))
    
    triangleListOutput = []
    
    pointIndex = 0
    while pointIndex < len(pointList):       
        currentPoint = pointList[pointIndex].copy()
        # move 'completed' triangles to the output
        i = 0
        while (i < len(triangleList)): 
            myTriangle = triangleList[i]
            if (triangle.isCircumCircleLeft(myTriangle, currentPoint)
            and not triangle.hasVertexInDomain(myTriangle)): 
                if triangle.isInsideDTM(myTriangle, header, dtm):
                    triangleListOutput.append(myTriangle)
                triangleList.pop(i)
            else: i += 1
        
        if Delaunay(triangleList, currentPoint):
            pointIndex += 1
        else:
            pointList.pop(pointIndex)
                
    # remove all the triangles belonging to the boundary 
    # and move the remaining triangles to the output
    while (len(triangleList) > 0):
        if not triangle.hasVertexInDomain(triangleList[0]):
            if triangle.isInsideDTM(triangleList[0], header, dtm):
                triangleListOutput.append(triangleList[0])
        triangleList.pop(0)
     
    return(pointList, triangleListOutput)

