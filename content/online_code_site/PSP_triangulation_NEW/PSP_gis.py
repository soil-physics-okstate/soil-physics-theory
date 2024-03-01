#PSP_gis.py
import numpy as np
from math import sqrt, fabs, floor

NODATA = -9999
NOLINK = -1

class Cheader():
    nrRows = 0
    nrCols = 0
    nrPoints = 0        
    xllCorner = NODATA        
    yllCorner = NODATA
    cellSize = 0
    flag = NODATA
    zMin = NODATA
    zMax = NODATA
    dz = NODATA
    magnify = NODATA
    
    
def computeMinMaxGrid(header, grid):
    isFirstValue = True
    nrValidPoints = 0
    for col in range(header.nrCols):
        for row in range(header.nrRows):
            z = grid[col, row]
            if (z != header.flag):
                nrValidPoints += 1
                if (isFirstValue):
                    header.zMin = z
                    header.zMax = z
                    isFirstValue = False
                else:
                    header.zMin = min(header.zMin, z)
                    header.zMax = max(header.zMax, z)
    
    header.nrPoints = nrValidPoints
    header.dz = max(1.0, header.zMax - header.zMin)                
    size = header.cellSize * sqrt(nrValidPoints)
    ratio = size / header.dz
    header.magnify = max(1., min(10., ratio / 4.))
    
    
def isOutOfGridColRow(header, col, row):
    if ((row < 0) or (row >= header.nrRows) 
    or (col < 0) or (col >= header.nrCols)): 
        return(True)
    else: 
        return(False)
    
def isOutOfGridXY(header, x, y):
    if ((x < header.xllCorner) or (y < header.yllCorner)
    or (x >= (header.xllCorner + header.cellSize * header.nrCols))
    or (y >= (header.yllCorner + header.cellSize * header.nrRows))):
        return(True)
    else: 
        return(False)
    
def getColRowFromXY(header, x, y):
    row = header.nrRows - floor((y - header.yllCorner) / header.cellSize)-1
    col = floor((x - header.xllCorner) / header.cellSize)
    return col, row;
    
#get value with check on boundary limits
def getValueFromColRow(header, grid, col, row):
    if (isOutOfGridColRow(header, col, row)): 
        return header.flag
    else: 
        return grid[col, row]
    
#get value with check on boundary limits         
def getValueFromXY(header, grid, x, y):
    if (isOutOfGridXY(header, x, y)): 
        return header.flag
    else:
        col, row = getColRowFromXY(header, x, y)
        return grid[col, row]
    
def getPointFromColRow(header, grid, col, row):
    x = header.xllCorner + (col + 0.5) * header.cellSize
    y = header.yllCorner + (header.nrRows - row - 0.5) * header.cellSize
    z = NODATA
    if not isOutOfGridColRow(header, col, row): 
        if grid[col, row] != header.flag:
            z = grid[col, row]
    return np.array([x,y,z])

def isTrue(header, booleanGrid, col, row):
    if isOutOfGridColRow(header, col, row): 
        return(False)
    else: 
        return booleanGrid[col][row]
  
def distance2D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    return sqrt(dx*dx + dy*dy) 

def distance3D(v1, v2):
    dx = fabs(v1[0] - v2[0])
    dy = fabs(v1[1] - v2[1])
    dz = fabs(v1[2] - v2[2])
    return sqrt(dx*dx + dy*dy + dz*dz)
