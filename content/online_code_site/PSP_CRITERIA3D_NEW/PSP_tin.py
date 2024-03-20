#PSP_tin.py
from math import fabs, sqrt
import numpy as np
from copy import copy
from PSP_public import NODATA, NOLINK
from IPython.testing.decorators import skip

class CheaderTin():     
    xMin = NODATA  
    xMax = NODATA
    yMin = NODATA
    yMax = NODATA
    zMin = NODATA
    zMax = NODATA 
    dz = NODATA     
    magnify = NODATA
    
class Ctriangle:
    def __init__(self, v = np.zeros((3, 3), float)):
        self.isBoundary = False
        self.boundarySlope = NODATA
        self.boundarySide = NODATA
        self.v = copy(v)
        if (not np.all(v == 0.0)): 
            self.centroid = (v[0]+v[1]+v[2])/3.0
            self.area = getArea2D(self.v)

#global structures
header = CheaderTin()
C3DTIN = []

def magnitude(v):
    return(np.sqrt(v.dot(v)))

def getArea(v):
    return 0.5 * magnitude(np.cross(v[1] - v[0], v[2] - v[0]))

# Area = 1/2 |x1(y2-y3) - x2(y1-y3) + x3(y1 - y2)| 
def getArea2D(v):
    x = v[:,0]
    y = v[:,1]
    return 0.5 * fabs(x[0]*(y[1]-y[2]) - x[1]*(y[0]-y[2]) + x[2]*(y[0]-y[1]))

def getHeader(triangleList):
    header = CheaderTin()
    header.xMin = triangleList[0].centroid[0]
    header.yMin = triangleList[0].centroid[1]
    header.zMin = triangleList[0].centroid[2]
    header.xMax = header.xMin
    header.yMax = header.yMin
    header.zMax = header.zMin
    
    for i in range(1, len(triangleList)):
        x = triangleList[i].centroid[0]
        y = triangleList[i].centroid[1]
        z = triangleList[i].centroid[2]
        header.xMin = min(header.xMin, x)
        header.yMin = min(header.yMin, y)
        header.zMin = min(header.zMin, z)
        header.xMax = max(header.xMax, x)
        header.yMax = max(header.yMax, y)
        header.zMax = max(header.zMax, z)
        
    dx = header.xMax - header.xMin
    dy = header.yMax - header.yMin
    header.dz = header.zMax - header.zMin
    ratio = sqrt(dx*dy) / header.dz
    header.magnify = max(1., min(10., ratio / 5.))
    return(header)
         
def distance2D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    return sqrt(dx*dx + dy*dy)

def distance3D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    dz = fabs(v1[2] - v2[2])
    return sqrt(dx*dx + dy*dy + dz*dz)


def boundaryProperties(TIN, index, neighbours):
    # compute nr of adj triangles 
    nrAdj = [0,0,0]
    for i in range(3):
        for j in range(3):
            if int(neighbours[j]) != NOLINK:
                for k in range(3):
                    if np.all (TIN[index].v[i] == TIN[int(neighbours[j])].v[k]):
                       nrAdj[i] += 1 
                       break
                   
    # internal and external vertices 
    internal = 0
    maxAdj = nrAdj[0]
    for i in np.arange(1, 3):
        if (nrAdj[i] > maxAdj) or (nrAdj[i] == maxAdj and TIN[index].v[i][2] > TIN[index].v[internal][2]):
            internal = i
            maxAdj = nrAdj[i]
    
    external = []
    for i in range(3):
        if (i != internal):
            external.append(TIN[index].v[i])
                       
    # compute boundary side
    TIN[index].boundarySide = distance3D(external[0], external[1])
    boundaryPoint = (external[0] + external[1]) * 0.5
    
    #compute slope
    dz = TIN[index].centroid[2] - boundaryPoint[2]
    dxy = distance2D(TIN[index].centroid, boundaryPoint)
    TIN[index].boundarySlope = dz/dxy


def getAdjacentVertices(t1, t2):
    isFirst = True
    for i in range(3):
        for j in range(3):
            if (t1[i] == t2[j]):
                if isFirst:
                    index1 = t1[i]
                    isFirst = False
                else:
                    index2= t1[i]
                    return (index1, index2)
    return NOLINK, NOLINK

def getAdjacentSide(i, j, vertexList, triangleList):
    triangle1 = triangleList[i]
    triangle2 = triangleList[j]
    index1, index2 = getAdjacentVertices(triangle1, triangle2)
    v1 = vertexList[int(index1)]
    v2 = vertexList[int(index2)]
    return distance2D(v1, v2)     

