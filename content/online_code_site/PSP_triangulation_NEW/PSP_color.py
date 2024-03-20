#PSP_gis.py
import numpy as np

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


def setColorScaleDtm():
    global colorScaleDtm
    keyColors = np.zeros((4,3),float)
    keyColors[0] = (0, 128, 0)      #green
    keyColors[1] = (255, 255, 0)    #yellow
    keyColors[2] = (128, 64, 0)     #brown
    keyColors[3] = (192, 192, 192)   #grey
    colorScaleDtm = setColorScale(256, keyColors)
    
def getDTMColor(z, header):
    zRelative = (z - header.zMin) / header.dz
    index = int(zRelative * (len(colorScaleDtm)-1))
    index = min(len(colorScaleDtm)-1, max(index, 0))
    return(colorScaleDtm[index])
