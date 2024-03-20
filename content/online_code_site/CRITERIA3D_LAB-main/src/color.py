# color.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

import numpy as np
from commonConst import EPSILON

colorScaleTIN = np.array([[], []], float)
colorRangeSE = np.array([[], []], float)
colorRangeRoot = np.array([[], []], float)
colorRangeSurfaceWater = np.array([[], []], float)


def setColorScale(nrLevels, keyColors):
    for i in range(len(keyColors)):
        for j in range(3):
            keyColors[i, j] /= 256.

    nrIntervals = len(keyColors) - 1
    step = int(max(nrLevels / nrIntervals, 1))
    myScale = np.zeros((nrLevels, 3), float)

    index = 0
    for i in range(nrIntervals):
        dRed = (keyColors[i + 1, 0] - keyColors[i, 0]) / step
        dGreen = (keyColors[i + 1, 1] - keyColors[i, 1]) / step
        dBlue = (keyColors[i + 1, 2] - keyColors[i, 2]) / step

        for j in range(step):
            index = step * i + j
            myScale[index, 0] = keyColors[i, 0] + (dRed * j)
            myScale[index, 1] = keyColors[i, 1] + (dGreen * j)
            myScale[index, 2] = keyColors[i, 2] + (dBlue * j)

    lastIndex = index
    if lastIndex < (nrLevels - 1):
        for i in range(lastIndex, nrLevels):
            myScale[i] = myScale[lastIndex]
    return myScale


def setColorScaleTIN():
    global colorScaleTIN
    keyColors = np.zeros((4, 3), float)
    keyColors[0] = (0, 128, 0)  # green
    keyColors[1] = (255, 255, 0)  # yellow
    keyColors[2] = (128, 64, 0)  # brown
    keyColors[3] = (192, 192, 192)  # grey
    colorScaleTIN = setColorScale(512, keyColors)


def setColorScaleDegreeOfSaturation():
    global colorRangeSE
    keyColors = np.zeros((4, 3), float)

    keyColors[0] = (255, 0, 0)  # red
    keyColors[1] = (255, 255, 0)  # yellow
    keyColors[2] = (0, 255, 0)  # green
    keyColors[3] = (0, 0, 255)  # blue
    colorRangeSE = setColorScale(1024, keyColors)


def setColorScaleRoot():
    global colorRangeRoot
    keyColors = np.zeros((4, 3), float)

    keyColors[0] = (128, 128, 128)  # gray
    keyColors[1] = (0, 255, 0)  # green
    keyColors[2] = (255, 255, 0)  # yellow
    keyColors[3] = (255, 0, 0)  # red

    colorRangeRoot = setColorScale(1024, keyColors)


def setColorScaleSurfaceWater():
    global colorRangeSurfaceWater
    keyColors = np.zeros((3, 3), float)
    keyColors[0] = (255, 255, 255)  # white
    keyColors[1] = (0, 255, 255)
    keyColors[2] = (0, 0, 255)  # blue
    colorRangeSurfaceWater = setColorScale(512, keyColors)


def setAllColorScale():
    setColorScaleTIN()
    setColorScaleDegreeOfSaturation()
    setColorScaleRoot()
    setColorScaleSurfaceWater()


def getTINColor(z, header):
    zRelative = 0.0 if (header.dz == 0.0) else (z - header.zMin) / header.dz
    index = int(zRelative * (len(colorScaleTIN) - 1))
    index = min(len(colorScaleTIN) - 1, max(index, 0))
    return colorScaleTIN[index]


def getSEColor(degreeSaturation, minimum, maximum):
    if degreeSaturation == 0:
        return 0.5, 0.5, 0.5
    percentage = (degreeSaturation - minimum) / (maximum - minimum)
    percentage = min(1.0, max(percentage, 0.0))
    index = int(percentage * (len(colorRangeSE) - 1))
    return colorRangeSE[index]


def getRootColor(rootDensity, minimum, maximum):
    value = (rootDensity - minimum) / (maximum - minimum)
    value = min(1.0, max(value, 0.0))
    index = int(value * (len(colorRangeSE) - 1))
    return colorRangeRoot[index]


# signPsi [m]
def getMatricPotentialColor(signPsi):
    signPsi *= 9.81     # [kPa]
    if signPsi > -0.1:
        return 0.5, 0, 1
    elif signPsi > -1.0:
        return 0, 0, 1
    elif signPsi > -5.0:
        return 0, 0.4, 1
    elif signPsi > -10:
        return 0, 0.8, 1
    elif signPsi > -20:
        return 0, 1, 0.8
    elif signPsi > -33:
        return 0, 1, 0
    elif signPsi > -50:
        return 0.8, 1, 0
    elif signPsi > -100:
        return 1, 0.75, 0
    elif signPsi > -300:
        return 1, 0.25, 0
    else:
        return 1, 0, 0


def getSurfaceWaterColor(waterHeight, maximum):
    percentage = waterHeight / maximum
    percentage = min(1.0, max(0.0, percentage))
    if percentage < EPSILON:
        return 0, 1, 0  # green
    else:
        # index = int(percentage * (len(colorRangeSE) - 1))
        # return colorRangeSE[index]
        return 0, 1 - percentage, percentage
