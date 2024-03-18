#PSP_travelTime.py
from __future__ import print_function, division
 
import math
import numpy as np

c = 299792458                      
NODATA = -9999
MAXDELTAINDEX = 8
SX = 0
DX = 1

class CLine:
    a = NODATA
    b = NODATA
    
class CPoint:
    x = NODATA
    y = NODATA
   
line1 = CLine()
line2 = CLine()
lastFlatLine = CLine() 

p0 = CPoint()
p1 = CPoint()
p2 = CPoint()
p3 = CPoint()

timeVector = []             
reflecCoeff = [] 
dy =[]
dy2 =[]

deltaSpace = 0             
deltaTime = 0               

def indexOfMaxVector(y, first, last):
    myMax = max(y[first:last])
    for i in range(first, last):
        if (y[i] == myMax): 
            return(i)

def indexOfMinVector(y, first, last):
    myMin = min(y[first:last])
    for i in range(first, last):
        if (y[i] == myMin): 
            return(i)
        
def getIndexOfFirstMinimum(y, first, last):
    threshold = -0.2
    indexMin = indexOfMinVector(y, first, last)
    index = indexMin
    valueMin = y[indexMin]
    indexMax = indexOfMaxVector(y, first, indexMin)
    
    while (valueMin < threshold):
        indexMin = indexOfMinVector(y, first, indexMax)
        valueMin = y[indexMin]
        if (valueMin < threshold):
            index = indexMin
            indexMax = indexOfMaxVector(y, first, indexMin)
    return(index)   
             
def avg(y, index1, index2):
    if (index2 < index1):return(NODATA)
    first = max(index1, 0)
    last = min(index2+1, len(y))
    nrValues = last - first
    return sum(y[first:last]) / nrValues

def normalizeVector(y):
    avgFirstValues = avg(y, 1, 6)
    return (y - avgFirstValues)
    
def WF_parameters(Vp, probeHandle, windowBegin, windowWidth, nrPoints):
    global deltaTime, deltaSpace, timeVector
    #abs. time [s] corresponding to the 1st point
    firstPointTime = 2. * windowBegin / (c*Vp)      
    deltaSpace = windowWidth /(nrPoints - 1)        
    deltaTime = 2. * deltaSpace / (c*Vp)             
    timeVector = np.zeros(nrPoints, float)            
    for i in range(nrPoints):
        timeVector[i] = firstPointTime + deltaTime * i
        
def runningAverage(y, nrPoints):   
    smooth = np.zeros(len(y), float) 
    for i in range(len(y)):
        smooth[i] = avg(y, i-nrPoints, i+nrPoints)
    return smooth / max(abs(smooth))

def firstDerivative5Points(y):
    d = np.zeros(len(y), float)
    for i in range(2):
        d[i] = 0.
    for i in range(2, len(y)-2):
        d[i] = (1./(12.)) * (y[i-2] - 8.*y[i-1] + 8.*y[i+1] - y[i+1])
    for i in range(len(y)-2, len(y)):
        d[i] = 0.      
    return d / max(abs(d))

# return a line structure with intercept (b) and slope (a) 
def weightedLinearRegression (x, y, index1, index2, versus):
    sumX = sumY = 0.
    sumX2 = sumXY = 0.
    
    if(index1 == index2):
        index1 -= 1
        index2 += 1
    
    #check index range
    if (index1 < 0):
        index1 = 0
    if (index2 >= len(y)):
        index2 = len(y)-1

    nrPoints = index2-index1+1
    if (versus == SX):
        for i in range(nrPoints-1, -1, -1): 
            for j in range (i+1):
                sumX += x[index1+i]
                sumY += y[index1+i]
                sumX2 += (x[index1+i]* x[index1+i])
                sumXY += x[index1+i] * y[index1+i]
    else:
        for i in range(nrPoints): 
            for j in range (i+1):
                sumX += x[index1+i]
                sumY += y[index1+i]
                sumX2 += (x[index1+i]* x[index1+i])
                sumXY += x[index1+i] * y[index1+i]
    
    n = (nrPoints*(nrPoints+1))/2
    line = CLine()
    line.a = (sumXY - sumX * sumY/n) / (sumX2 - sumX * sumX/n)
    line.b = (sumY - line.a * sumX)/n
    return(line)


#backward function 
def checkZeroValue(y, indexMaxY):
    index = indexMaxY
    while ((y[index] > 0.01) and (index > 0)):
        index -= 1
    if ((index == 0) and (y[index] > 0.01)):
        return(NODATA)
    else:
        if (abs(y[index]) < abs(y[index+1])):
            return (index)
        else:
            return (index+1)

def lineIntersection(line1, line2):
    myPoint = CPoint()
    if (line1.a != line2.a):
        myPoint.x = (line2.b - line1.b) / (line1.a - line2.a)
        myPoint.y = myPoint.x * line1.a + line1.b
    else:
        myPoint.x = NODATA
        myPoint.y = NODATA
    
    return(myPoint)

def computeTravelTime(probeHandle, permittivity, Vp):
    global dy, dy2, line1, line2, lastFlatLine
    global indexFlatLine, indexRegr1, indexRegr2, indexRegr3
    global p0, p1, p2, p3
    
    dy = firstDerivative5Points(reflecCoeff)
    dy2 = firstDerivative5Points(dy)
    dy = runningAverage(dy, 3)
    dy2 = runningAverage(dy2, 3)
    
    indexP1 = getIndexOfFirstMinimum(dy2, 0, len(dy2))
    indexP0 = indexOfMaxVector(dy2, 0, indexP1)
    # p0 or p1 not found
    if indexP0 == 0 or indexP1 == 0:
        return False
    
    p0.x = indexP0 * deltaTime
    p0.y = reflecCoeff[indexP0]
    
    p1.x = indexP1 * deltaTime
    p1.y = reflecCoeff[indexP1]

    #search second reflection
    indexMaxDerivative = indexOfMaxVector(dy, indexP1, len(dy))
    indexZeroDerivative = checkZeroValue(dy, indexMaxDerivative)
    delta = min((indexMaxDerivative - indexZeroDerivative), MAXDELTAINDEX)
    indexRegr2 = indexZeroDerivative - delta
    indexRegr3 = indexZeroDerivative + delta
    
    line1 = weightedLinearRegression(timeVector, reflecCoeff, 
                        indexRegr2 - delta, indexRegr2 + delta, DX)
    
    line2 = weightedLinearRegression(timeVector, reflecCoeff, 
                        indexRegr3 - delta, indexRegr3 + delta, SX)
    
    p2 = lineIntersection(line1, line2)
    index = int(p2.x / deltaTime)
    p2.y = reflecCoeff[index]
    
    #compute p3 (asymptote)
    last = len(reflecCoeff)
    p3.x = (last-1) * deltaTime
    p3.y = avg(reflecCoeff, last-MAXDELTAINDEX, last-1) 
    lastFlatLine = CLine()
    lastFlatLine.a = 0;
    lastFlatLine.b = p3.y   
    return True
