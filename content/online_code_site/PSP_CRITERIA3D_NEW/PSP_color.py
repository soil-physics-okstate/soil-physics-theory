#PSP_color.py
import numpy as np

colorScaleTIN = np.array([[],[]], float)
colorRangeSE = np.array([[],[]], float)
colorRangeSurfaceWater = np.array([[],[]], float)

def setColorScale(nrLevels, keyColors):
    for i in range (len(keyColors)):
        for j in range (3):
            keyColors[i,j] /= 256.
            
    nrIntervals = len(keyColors)-1
    step = int(max(nrLevels / nrIntervals, 1))
    myScale = np.zeros((nrLevels, 3), float)
    
    for i in range (nrIntervals):
        dRed = (keyColors[i+1,0] - keyColors[i,0]) / step
        dGreen = (keyColors[i+1,1] - keyColors[i,1]) / step
        dBlue = (keyColors[i+1,2] - keyColors[i,2]) / step
   
        for j in range (step):
            index = step * i + j
            myScale[index, 0] = keyColors[i,0] + (dRed * j)
            myScale[index, 1] = keyColors[i,1] + (dGreen * j)
            myScale[index, 2] = keyColors[i,2] + (dBlue * j)
    
    lastIndex = index
    if (lastIndex < (nrLevels-1)):
        for i in range(lastIndex, nrLevels):
            myScale[i] = myScale[lastIndex]
    return (myScale)

def setColorScaleTIN():
    global colorScaleTIN
    keyColors = np.zeros((4,3),float)
    keyColors[0] = (0, 128, 0)      #green
    keyColors[1] = (255, 255, 0)    #yellow
    keyColors[2] = (128, 64, 0)     #brown
    keyColors[3] = (192, 192, 192)   #grey
    colorScaleTIN = setColorScale(512, keyColors)

def setColorScaleDegreeOfSaturation():
    global colorRangeSE
    keyColors = np.zeros((3,3),float)
    
    keyColors[0] = (255, 0, 0)       # red
    keyColors[1] = (255, 255, 0)     # yellow
    keyColors[2] = (0, 0, 255)       # blue
    colorRangeSE = setColorScale(1024, keyColors)
    
def setColorScaleSurfaceWater():
    global colorRangeSurfaceWater
    keyColors = np.zeros((3,3),float)
    keyColors[0] = (255, 255, 255)   #white
    keyColors[1] = (0, 255, 255)     
    keyColors[2] = (0, 0, 255)       #blue
    colorRangeSurfaceWater = setColorScale(512, keyColors)
    
def setAllColorScale():
    setColorScaleTIN()
    setColorScaleDegreeOfSaturation()
    setColorScaleSurfaceWater()
    
def getTINColor(z, header):
    zRelative = (z - header.zMin) / header.dz
    index = int(zRelative * (len(colorScaleTIN)-1))
    index = min(len(colorScaleTIN)-1, max(index, 0))
    return(colorScaleTIN[index])

def getSEColor(degreeSaturation, minimum, maximum):
    perc = (degreeSaturation-minimum) / (maximum-minimum)
    perc = min(1.0, max(0.0, perc))
    index = int(perc * (len(colorRangeSE)-1))
    return(colorRangeSE[index])

def getSurfaceWaterColor(waterHeight, maximum):
    perc = waterHeight / maximum
    perc = min(1.0, max(0.0, perc))
    index = int(perc*(len(colorRangeSurfaceWater)-1))
    return colorRangeSurfaceWater[index]

