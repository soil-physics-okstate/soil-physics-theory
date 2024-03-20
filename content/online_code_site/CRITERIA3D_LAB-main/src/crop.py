# crop.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------


from dataStructures import *
import math
import soil
import rectangularMesh
import numpy as np

MAX_EVAPORATION_DEPTH = 0.15  # [m]


class CCrop:
    laiMonth = [float]  # [m2 m-2]
    rootDepthZero = NODATA  # [m]
    rootDepthMax = NODATA  # [m]
    rootWidth = NODATA  # [m]
    rootXDeformation = NODATA  # [-]
    rootZDeformation = NODATA  # [-]
    kcMax = NODATA  # [-]
    fRAW = NODATA  # [-]
    currentLAI = NODATA
    currentRootDepth = NODATA
    currentRootLength = NODATA

    def setMaxRootDepth(self):
        self.currentRootDepth = self.rootDepthMax
        self.currentRootLength = self.currentRootDepth - self.rootDepthZero

    def setCurrentLAI(self, currentDate):
        m = currentDate.month
        lai0 = float(self.laiMonth[m - 1])
        if m < 12:
            lai1 = float(self.laiMonth[m])
        else:
            lai1 = float(self.laiMonth[0])
        # TO DO improve for month length
        self.currentLAI = lai0 + (lai1 - lai0) * (currentDate.day - 1) / 30


# global variables
currentCrop = CCrop()
rootDensity = []
k_root = np.array([], np.float64)
global maxRootFactor, x, y


# initialize root density and root factor to zero
def setNoCrop():
    global k_root, rootDensity, maxRootFactor

    k_root = np.zeros(C3DStructure.nrRectangles)
    rootDensity.clear()
    for i in range(C3DStructure.nrRectangles):
        rootDensity.append(np.zeros(C3DStructure.nrLayers, np.float64))
    maxRootFactor = 1


def initializeCrop():
    global k_root, rootDensity, maxRootFactor
    global x, y

    # check maximum root depth
    if currentCrop.rootDepthMax > C3DStructure.gridDepth:
        print("WARNING!\ncurrentCrop.rootDepthMax > gridDepth")
        print("currentCrop.rootDepthMax will be set = " + str(C3DStructure.gridDepth))
        currentCrop.rootDepthMax = C3DStructure.gridDepth

    currentCrop.setMaxRootDepth()

    # initialize root factor
    k_root = np.zeros(C3DStructure.nrRectangles)

    # distance from plant row (x axis)
    max_distance = currentCrop.rootWidth * 0.5
    a = y[1] - y[0]
    b = x[0] - x[1]
    c = y[0] * (x[1] - x[0]) - x[0] * (y[1] - y[0])
    denominator = math.sqrt(a * a + b * b)

    for i in range(C3DStructure.nrRectangles):
        [cx, cy, cz] = rectangularMesh.C3DRM[i].centroid
        line_distance = math.fabs(a * cx + b * cy + c) / denominator

        if line_distance < max_distance or math.fabs(line_distance - max_distance) < EPSILON:
            factor = 1.0 - line_distance / (max_distance * 0.5)
            k_root[i] = 1.0 + factor * currentCrop.rootXDeformation
        else:
            k_root[i] = 0.0

    # set root density
    rootDensity.clear()
    maxRootFactor = 0
    for i in range(C3DStructure.nrRectangles):
        rootDensity.append(computeRootDensity(currentCrop, C3DStructure.nrLayers, k_root[i]))
        # update max root factor
        for layer in range(1, C3DStructure.nrLayers):
            root_factor = k_root[i] * rootDensity[i][layer] / soil.thickness[layer]
            maxRootFactor = max(root_factor, maxRootFactor)


# fraction of intercepted photosynthetically active radiation [-]
def fPARi(currentLAI):
    ke = 0.6  # light extinction coefficient [-]
    if (currentLAI == NODATA) or (currentLAI <= 0):
        return 0.
    else:
        # Beer–Lambert Equation
        return 1. - math.exp(-ke * currentLAI)


def getMaxEvaporation(currentLAI, ET0):
    return ET0 * 0.66 * (1. - fPARi(currentLAI))


def getMaxTranspiration(currentLAI, kcMax, ET0):
    if (currentLAI == NODATA) or (currentLAI <= 0):
        return 0.
    else:
        fPAR = fPARi(currentLAI)
        TC = 1 + (kcMax - 1) * fPAR
        return ET0 * fPAR * TC


def cardioidDistribution(deformationFactor, nrLayersWithRoot):
    # initialize
    halfCircle = np.zeros(nrLayersWithRoot)
    cardioid = np.zeros(nrLayersWithRoot * 2)

    for i in range(nrLayersWithRoot):
        sinAlfa = 1.0 - float(i + 1.0) / float(nrLayersWithRoot)
        cosAlfa = max(math.sqrt(1.0 - math.pow(sinAlfa, 2.0)), EPSILON)
        alfa = math.atan(sinAlfa / cosAlfa)
        halfCircle[i] = ((math.pi / 2.0) - alfa - sinAlfa * cosAlfa) / math.pi

    lastLayer = nrLayersWithRoot * 2 - 1
    cardioid[0] = halfCircle[0]
    cardioid[lastLayer] = cardioid[0]
    for i in range(1, nrLayersWithRoot):
        cardioid[i] = halfCircle[i] - halfCircle[i - 1]
        cardioid[lastLayer - i] = cardioid[i]

    # cardioid deformation
    LiMin = -math.log(0.2) / float(nrLayersWithRoot)
    LiMax = -math.log(0.05) / float(nrLayersWithRoot)
    k = LiMin + (LiMax - LiMin) * (deformationFactor - 1.0)

    rootDensitySum = 0
    for i in range(nrLayersWithRoot * 2):
        cardioid[i] *= math.exp(-k * (i + 0.5))
        rootDensitySum += cardioid[i]

    # normalize
    for i in range(nrLayersWithRoot * 2):
        cardioid[i] /= rootDensitySum

    # assign layer density
    layerDensity = np.zeros(nrLayersWithRoot, np.float64)
    for i in range(nrLayersWithRoot):
        layerDensity[i] = cardioid[2 * i] + cardioid[2 * i + 1]

    return layerDensity


def computeRootDensity(crop, nrLayers, rootFactor):
    # initialize
    myRootDensity = np.zeros(nrLayers, np.float64)
    if crop.currentRootLength <= 0 or rootFactor == 0:
        return myRootDensity

    rootLength = crop.currentRootLength * min(1.0, math.sqrt(rootFactor))
    rootZero = crop.rootDepthZero
    if rootLength < 0.001:
        return myRootDensity

    # smallest unit of computation (1 mm)
    atoms = np.zeros(nrLayers, np.int16)
    for i in range(nrLayers):
        atoms[i] = int(round(soil.thickness[i] * 1000))

    nrUnrootedAtoms = int(round(rootZero * 1000))
    nrRootedAtoms = int(round(rootLength * 1000))
    densityAtoms = cardioidDistribution(crop.rootZDeformation, nrRootedAtoms)

    # assign root density
    counter = 0
    for layer in range(nrLayers):
        for i in range(atoms[layer]):
            if (counter >= nrUnrootedAtoms) and (counter < nrUnrootedAtoms + nrRootedAtoms):
                myRootDensity[layer] += densityAtoms[counter - nrUnrootedAtoms]
            counter += 1

    # check (rootDensitySum == 1)
    rootDensitySum = 0.0
    for i in range(nrLayers):
        rootDensitySum += myRootDensity[i]
    if abs(rootDensitySum - 1.0) > EPSILON:
        print("WARNING! Sum of root density:", rootDensitySum)

    return myRootDensity


# assign hourly transpiration
def setTranspiration(surfaceIndex, myRootDensity, maxTranspiration):
    if maxTranspiration < EPSILON:
        return 0.0

    # Initialize
    rootDensityWithoutStress = 0.0  # [-]
    actualTranspiration = 0.0  # [mm]

    nrLayers = len(myRootDensity)
    isLayerStressed = np.zeros(nrLayers, dtype=bool)
    layerTranspiration = np.zeros(nrLayers, np.float64)  # [mm]

    for layer in range(nrLayers):
        if myRootDensity[layer] > 0:
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            theta = soil.getVolumetricWaterContent(i)
            horizon = soil.horizons[C3DCells[i].horizonIndex]
            # water surplus
            if theta > horizon.FC:
                fraction = 1.0 - (theta - horizon.FC) / (horizon.thetaS - horizon.FC)
                fraction = fraction ** 3
                layerTranspiration[layer] = maxTranspiration * fraction * myRootDensity[layer]
                isLayerStressed[layer] = True
            else:
                # water scarcity
                wsThreshold = horizon.FC - currentCrop.fRAW * (horizon.FC - horizon.WP)
                if theta < wsThreshold:
                    if theta <= horizon.WP:
                        layerTranspiration[layer] = 0.0
                    else:
                        layerTranspiration[layer] = maxTranspiration * myRootDensity[layer] * (
                                (theta - horizon.WP) / (wsThreshold - horizon.WP))
                    isLayerStressed[layer] = True
                else:
                    # normal conditions
                    layerTranspiration[layer] = maxTranspiration * myRootDensity[layer]

                    # check stress
                    theta_mm = theta * soil.thickness[layer] * 1000.
                    WSThreshold_mm = wsThreshold * soil.thickness[layer] * 1000.

                    if (theta_mm - layerTranspiration[layer]) > WSThreshold_mm:
                        isLayerStressed[layer] = False
                        rootDensityWithoutStress += myRootDensity[layer]
                    else:
                        isLayerStressed[layer] = True

            actualTranspiration += layerTranspiration[layer]

    # WATER STRESS [-]
    waterStress = 1. - (actualTranspiration / maxTranspiration)

    # Hydraulic redistribution
    # the movement of water from moist to dry soil through plant roots
    # TODO add numerical process
    if waterStress > EPSILON and rootDensityWithoutStress > EPSILON:
        redistribution = min(waterStress, rootDensityWithoutStress) * maxTranspiration

        # redistribution acts on not stressed roots
        for layer in range(nrLayers):
            if (myRootDensity[layer] > 0) and (not isLayerStressed[layer]):
                addTranspiration = redistribution * (myRootDensity[layer] / rootDensityWithoutStress)
                layerTranspiration[layer] += addTranspiration
                actualTranspiration += addTranspiration

    # Assign transpiration flux [m3 s-1]
    for layer in range(nrLayers):
        if layerTranspiration[layer] > 0:
            i = surfaceIndex + C3DStructure.nrRectangles * layer
            rate = (layerTranspiration[layer] * 0.001) / 3600.0  # [m s-1]
            C3DCells[i].sinkSource -= rate * C3DCells[i].area  # [m3 s-1]

    return actualTranspiration


# assign hourly soil evaporation
def setEvaporation(surfaceIndex, maxEvaporation):
    # TODO: enable surface evaporation - numerical problem
    """
    surfaceWater = (C3DCells[surfaceIndex].H - C3DCells[surfaceIndex].z)        # [m]
    surfaceEvaporation = min(maxEvaporation, surfaceWater * 1000.0)             # [mm]
    rate = (surfaceEvaporation * 0.001) / 3600.0                                # [m s-1]
    C3DCells[surfaceIndex].sinkSource -= rate * C3DCells[surfaceIndex].area     # [m3 s-1]
    """
    surfaceEvaporation = 0
    actualEvaporation = surfaceEvaporation
    residualEvaporation = maxEvaporation - surfaceEvaporation

    if residualEvaporation < EPSILON:
        return actualEvaporation

    # soil evaporation
    lastIndex = 0
    while (lastIndex < len(soil.depth)) and (soil.depth[lastIndex] <= MAX_EVAPORATION_DEPTH):
        lastIndex += 1

    nrEvapLayers = lastIndex

    # depth coefficient: 1 at first layer, ~0.1 at MAX_EVAPORATION_DEPTH
    evapCoefficient = np.zeros(nrEvapLayers, np.float64)
    sumCoefficient = 0
    for i in range(1, nrEvapLayers):
        depthCoefficient = max((soil.depth[i] - soil.depth[1]) / (MAX_EVAPORATION_DEPTH - soil.depth[1]), 0)
        evapCoefficient[i] = math.exp(-depthCoefficient * math.e)
        sumCoefficient += (evapCoefficient[i] * soil.thickness[i])

    isWaterSupply = True
    while (residualEvaporation > EPSILON) and (isWaterSupply is True):
        isWaterSupply = False
        sumEvaporation = 0.0

        for layer in range(1, nrEvapLayers):
            index = surfaceIndex + C3DStructure.nrRectangles * layer
            horizon = soil.horizons[C3DCells[index].horizonIndex]
            theta = soil.getVolumetricWaterContent(index)  # [m3 m-3]
            half_FC = horizon.HYGR + (horizon.FC - horizon.HYGR) * 0.5
            evaporationThreshold = half_FC - evapCoefficient[layer] * (half_FC - horizon.HYGR)  # [m3 m-3]
            evaporation = residualEvaporation * ((evapCoefficient[layer]
                                                  * soil.thickness[layer]) / sumCoefficient)  # [mm]
            if theta < evaporationThreshold:
                evaporation = 0.0
            else:
                availableWC = (theta - evaporationThreshold) * soil.thickness[layer] * 1000.  # [mm]
                if availableWC <= evaporation:
                    evaporation = availableWC
                else:
                    isWaterSupply = True

            sumEvaporation += evaporation

            rate = (evaporation * 0.001) / 3600.  # [m s-1]
            C3DCells[index].sinkSource -= rate * C3DCells[index].area  # [m3 s-1]

        residualEvaporation -= sumEvaporation
        actualEvaporation += sumEvaporation

    return actualEvaporation


def setEvapotranspiration(currentDate, ET0):
    if C3DParameters.computeTranspiration:
        currentCrop.setCurrentLAI(currentDate)
        for i in range(C3DStructure.nrRectangles):
            maxTranspiration = getMaxTranspiration(currentCrop.currentLAI, currentCrop.kcMax, ET0) * k_root[i]
            setTranspiration(i, rootDensity[i], maxTranspiration)
    else:
        currentCrop.currentLAI = 0

    if C3DParameters.computeEvaporation:
        maxEvaporation = getMaxEvaporation(currentCrop.currentLAI, ET0)
        for i in range(C3DStructure.nrRectangles):
            setEvaporation(i, maxEvaporation)
