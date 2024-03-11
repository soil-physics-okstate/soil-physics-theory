# rectangularMesh.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from math import fabs, sqrt
import numpy as np
import soil
from copy import copy
from dataStructures import *
from enum import Enum


class Neighbours(Enum):
    Right = 1
    Up = 2
    Left = 3
    Down = 4


class CRectangularMeshHeader:
    xMin = NODATA
    xMax = NODATA
    yMin = NODATA
    yMax = NODATA
    zMin = NODATA
    zMax = NODATA
    dz = NODATA
    magnify = NODATA


# global structures
header = CRectangularMeshHeader()
C3DRM = []


class CRectangle:
    def __init__(self, index, v=np.zeros((C3DStructure.nrVerticesPerRectangle, C3DStructure.nrDimensions), float)):
        self.index = index
        self.v = copy(v)
        self.centroid = getCentroid3D(self.v)
        self.area = getArea2D(self.v)

        # initialize
        self.isBoundary = False
        self.boundarySide = NODATA
        self.boundarySlope = NODATA

        self.neighbours, self.boundarySides, self.boundarySlopes = boundaryProperties(self.index, self.v, self.centroid)

        if NOLINK in self.neighbours:
            isFirst = True
            index = NOLINK
            for i in range(len(self.boundarySlopes)):
                if self.boundarySlopes[i] != NODATA and self.boundarySlopes[i] > 0:
                    if isFirst:
                        index = i
                        isFirst = False
                    else:
                        if self.boundarySlopes[i] > self.boundarySlopes[index]:
                            index = i
            if index != NOLINK:
                self.isBoundary = True
                self.boundarySide = self.boundarySides[index]
                self.boundarySlope = self.boundarySlopes[index]


def rectangularMeshCreation():
    global header
    C3DRM.clear()
    C3DStructure.totalArea = 0
    dx = C3DStructure.gridWidth * 0.5
    dy = C3DStructure.gridHeight * 0.5

    index = 0
    for y in np.arange(-dy, dy, C3DStructure.cellSize):
        for x in np.arange(-dx, dx, C3DStructure.cellSize):
            if (x < dx) and (y < dy):
                rectangle = CRectangle(index, getRectangleVertices(x, y))
                C3DRM.append(copy(rectangle))
                C3DStructure.totalArea += rectangle.area
                index += 1
    header = getHeader(C3DRM)


def getRectangleVertices(x, y):
    v = np.zeros((C3DStructure.nrVerticesPerRectangle, C3DStructure.nrDimensions), float)
    dzy = C3DStructure.slopeY * y
    dzy2 = C3DStructure.slopeY * (y + C3DStructure.cellSize)
    dzPlant = C3DStructure.plantSlope * C3DStructure.plantSlopeWidth
    dzPlant1 = C3DStructure.plantSlope * min(C3DStructure.plantSlopeWidth, fabs(x))
    dzPlant2 = C3DStructure.plantSlope * min(C3DStructure.plantSlopeWidth, fabs(x + C3DStructure.cellSize))
    v[0] = [x, y, C3DStructure.elevation + dzPlant - dzPlant1 - dzy]
    v[1] = [x + C3DStructure.cellSize, y, C3DStructure.elevation + dzPlant - dzPlant2 - dzy]
    v[2] = [x + C3DStructure.cellSize, y + C3DStructure.cellSize, C3DStructure.elevation + dzPlant - dzPlant2 - dzy2]
    v[3] = [x, y + C3DStructure.cellSize, C3DStructure.elevation + dzPlant - dzPlant1 - dzy2]
    return v


def boundaryProperties(index, vertices, centroid):
    neighbours = []
    boundarySides = []
    boundarySlopes = []
    greater_than = lambda a, b: a >= b
    less_than = lambda a, b: a <= b
    plus = lambda a, b: a + b
    minus = lambda a, b: a - b
    neighbours_map = {
        Neighbours.Right: {
            "index": index + 1,
            "limit": C3DStructure.gridWidth*0.5,
            "comparator": less_than,
            "vertices": [vertices[1], vertices[2]],
            "boundary_coordinate": 0,
            "operator": plus
        },
        Neighbours.Up: {
            "index": index + C3DStructure.nrRectanglesInXAxis,
            "limit": C3DStructure.gridHeight*0.5,
            "comparator": less_than,
            "vertices": [vertices[2], vertices[3]],
            "boundary_coordinate": 1,
            "operator": plus
        },
        Neighbours.Left: {
            "index": index - 1,
            "limit": -C3DStructure.gridWidth*0.5,
            "comparator": greater_than,
            "vertices": [vertices[3], vertices[0]],
            "boundary_coordinate": 0,
            "operator": minus
        },
        Neighbours.Down: {
            "index": index - C3DStructure.nrRectanglesInXAxis,
            "limit": -C3DStructure.gridHeight*0.5,
            "comparator": greater_than,
            "vertices": [vertices[0], vertices[1]],
            "boundary_coordinate": 1,
            "operator": minus
        },
    }
    for key, value in neighbours_map.items():
        coordinate_to_check = neighbours_map[key]["boundary_coordinate"]
        limit_to_check = neighbours_map[key]["limit"]
        move_coordinate = neighbours_map[key]["operator"]
        check_limit = neighbours_map[key]["comparator"]

        if check_limit(move_coordinate(
                centroid[coordinate_to_check], C3DStructure.cellSize),
                limit_to_check):
            neighbours.append(neighbours_map[key]["index"])
            boundarySides.append(NODATA)
            boundarySlopes.append(NODATA)
        else:
            neighbours.append(NOLINK)

            # compute boundary side
            boundarySides.append(distance3D(neighbours_map[key]["vertices"][0], neighbours_map[key]["vertices"][1]))
            boundaryPoint = (neighbours_map[key]["vertices"][0] + neighbours_map[key]["vertices"][1]) * 0.5

            # compute slope
            dz = centroid[2] - boundaryPoint[2]
            dxy = distance2D(centroid, boundaryPoint)
            if abs(dz / dxy) < EPSILON:
                boundarySlopes.append(0.0)
            else:
                boundarySlopes.append(dz / dxy)
    return neighbours, boundarySides, boundarySlopes


def getCentroid3D(vertices):
    return np.array(
        [(vertices[0][0] + vertices[1][0]) / 2,
         (vertices[0][1] + vertices[2][1]) / 2,
         (vertices[0][2] + vertices[1][2] + vertices[2][2] + vertices[3][2]) / 4])


def getArea2D(vertices):
    x_side = vertices[1][0] - vertices[0][0]
    y_side = vertices[2][1] - vertices[0][1]
    return x_side * y_side


def getHeader(rectangleList):
    myHeader = CRectangularMeshHeader()
    myHeader.xMin = rectangleList[0].centroid[0]
    myHeader.yMin = rectangleList[0].centroid[1]
    myHeader.zMin = rectangleList[0].centroid[2]
    myHeader.xMax = myHeader.xMin
    myHeader.yMax = myHeader.yMin
    myHeader.zMax = myHeader.zMin

    for i in range(1, len(rectangleList)):
        x = rectangleList[i].centroid[0]
        y = rectangleList[i].centroid[1]
        z = rectangleList[i].centroid[2]
        myHeader.xMin = min(myHeader.xMin, x)
        myHeader.yMin = min(myHeader.yMin, y)
        myHeader.zMin = min(myHeader.zMin, z)
        myHeader.xMax = max(myHeader.xMax, x)
        myHeader.yMax = max(myHeader.yMax, y)
        myHeader.zMax = max(myHeader.zMax, z)

    halfStep = C3DStructure.cellSize * 0.5
    myHeader.xMin -= halfStep
    myHeader.xMax += halfStep
    myHeader.yMin -= halfStep
    myHeader.yMax += halfStep

    dx = myHeader.xMax - myHeader.xMin
    dy = myHeader.yMax - myHeader.yMin
    myHeader.magnify = max(1., 1. / sqrt(max(dx, dy)))
    return myHeader


def distance2D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    return sqrt(dx * dx + dy * dy)


def distance3D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    dz = fabs(v1[2] - v2[2])
    return sqrt(dx * dx + dy * dy + dz * dz)


def getAdjacentVertices(rectangle1, rectangle2):
    commonVertices = []
    for i in rectangle1.v:
        for j in rectangle2.v:
            if np.allclose(i, j, rtol=1e-05):
                commonVertices.append(i)
    if len(commonVertices) > 2:
        raise Exception('Too many common vertices')
    return commonVertices[0], commonVertices[1]


def getAdjacentSide(i, j):
    rectangle1 = C3DRM[i]
    rectangle2 = C3DRM[j]
    v1, v2 = getAdjacentVertices(rectangle1, rectangle2)
    return distance2D(v1, v2)


def isInsideRectangle(x, y, rectangle):
    if rectangle.v[0][0] <= x < rectangle.v[1][0] \
       and rectangle.v[1][1] <= y < rectangle.v[2][1]:
        return True
    else:
        return False


def getSurfaceIndex(x, y):
    # search index
    for index in range(C3DStructure.nrRectangles):
        if isInsideRectangle(x, y, C3DRM[index]):
            return index
    return NODATA


def getCellIndex(x, y, depth):
    surfaceIndex = getSurfaceIndex(x, y)
    if surfaceIndex != NODATA:
        for layer in range(C3DStructure.nrLayers):
            top = soil.depth[layer] - (soil.thickness[layer] * 0.5)
            bottom = soil.depth[layer] + (soil.thickness[layer] * 0.5)
            if top < depth <= bottom:
                index = surfaceIndex + C3DStructure.nrRectangles * layer
                return index
    return NODATA


def getXYDepth(index):
    x = C3DCells[index].x
    y = C3DCells[index].y
    surfaceIndex = getSurfaceIndex(x, y)
    depth = C3DCells[surfaceIndex].z - C3DCells[index].z
    return x, y, depth
