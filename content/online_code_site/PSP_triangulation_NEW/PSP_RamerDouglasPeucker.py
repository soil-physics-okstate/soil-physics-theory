#PSP_RamerDouglasPeucker.py
from __future__ import print_function, division
from PSP_gis import *

# return first boundary point, if exist
# scan grid from [0,0]
def getFirstPoint(header, boundary):
    for col in range(0, header.nrCols):
        for row in range(0, header.nrRows):
            if (boundary[col, row] == True):
                return(True, col, row)
    return(False, NODATA, NODATA)

# return first boundary point, if exist 
# in the neighbours of [row,col], anti-clockwise order     
def getNearestPoint(header, boundary, col, row):
    if isTrue(header, boundary, col, row+1):
        return(True, col, row+1)
    elif isTrue(header, boundary, col+1, row+1):
        return(True, col+1, row+1)
    elif isTrue(header, boundary, col+1, row):
        return(True, col+1, row)
    elif isTrue(header, boundary, col+1, row-1):
        return(True, col+1, row-1)
    elif isTrue(header, boundary, col, row-1):
        return(True, col, row-1)
    elif isTrue(header, boundary, col-1, row-1):
        return(True, col-1, row-1)
    elif isTrue(header, boundary, col-1, row):
        return(True, col-1, row)
    elif isTrue(header, boundary, col-1, row+1):
        return(True, col-1, row+1)
    return(False, NODATA, NODATA) 

def triangleAreaErone(v1, v2, v3):
    a = distance3D(v2, v3)
    b = distance3D(v1, v3)
    c = distance3D(v1, v2)
    #semiperimeter
    sp = (a + b + c) * 0.5
    squaredArea = sp*(sp-a)*(sp-b)*(sp-c)
    if (squaredArea < 0.0):
        return(0.0)
    else: 
        return sqrt(squaredArea)

def perpendicularDistance(p, v0, v1, base):
    area = triangleAreaErone(p,v0,v1)
    return (2.0 * area) / base

def firstBoundaryPoint(header, boundary):
    for col in range(0, header.nrCols):
        for row in range(0, header.nrRows):
            if (boundary[col, row] == True):
                return(True, col, row)
    return(False, NODATA, NODATA)


# reduce boundary points with Ramer-Douglas-Peucker algorithm
def RamerDouglasPeucker(pointList, threshold, maxSide):
    if len(pointList) < 3: return pointList
    dmax = 0
    index = 0
    end = len(pointList)-1
    
    base = distance3D(pointList[0], pointList[end])
    if (base > maxSide):
        index = int(end / 2.)
        dmax = base
    else:
        # Find the point with the maximum distance
        for i in range (1, end):
            d = perpendicularDistance(pointList[i], pointList[0], pointList[end], base) 
            if (d > dmax):
                index = i
                dmax = d
                
    # if distance > threshold cut the array 
    if (dmax > threshold):
        list1 = RamerDouglasPeucker(pointList[0:index+1], threshold, maxSide)
        list2 = RamerDouglasPeucker(pointList[index:len(pointList)], threshold, maxSide)
        
        resultList = []
        for j in range(len(list1)):   
            resultList.append(list1[j]) 
        for j in range(1, len(list2)):     
            resultList.append(list2[j])
    else:
        resultList = []
        resultList.append(pointList[0])
        resultList.append(pointList[end])
 
    return resultList 


def getDtmBoundary(header, dtm, threshold, maxSide):
    pointList = []
    boundary = np.zeros((header.nrCols, header.nrRows), bool)
    for col in range(header.nrCols):
        for row in range(header.nrRows):
            boundary[col, row] = False
            if dtm[col, row] != header.flag:
                if ((row == 0) or (col == 0) 
                or (row == (header.nrRows-1)) or (col == (header.nrCols-1))):
                    boundary[col, row] = True
                if row > 0 and dtm[col, row-1] == header.flag:
                    boundary[col, row] = True
                elif col > 0 and dtm[col-1, row] == header.flag:
                    boundary[col, row] = True
                elif row < (header.nrRows-1) and dtm[col, row+1] == header.flag:
                    boundary[col, row] = True
                elif col < (header.nrCols-1) and dtm[col+1, row] == header.flag:
                    boundary[col, row] = True
                    
    print ("reduce boundary points (Ramer-Douglas-Peucker algorithm)...")
        
    isFirstBoundary, col, row = getFirstPoint(header, boundary)
    while (isFirstBoundary):                    
        boundaryList = []
        isBoundary = True
        while (isBoundary):
            myPoint = getPointFromColRow(header, dtm, col, row)
            boundaryList.append(myPoint)
            boundary[col, row] = False  
            isBoundary, col, row = getNearestPoint(header, boundary, col, row)
       
        reducedList = RamerDouglasPeucker(boundaryList, threshold, maxSide) 
        
        for i in range(len(reducedList)): 
            pointList.append(reducedList[i])
             
        isFirstBoundary, col, row = getFirstPoint(header, boundary) 
        
    return(pointList)
