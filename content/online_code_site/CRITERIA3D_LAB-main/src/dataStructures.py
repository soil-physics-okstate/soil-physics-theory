# dataStructures.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from commonConst import *


class C3DStructure:
    nrDimensions = 3
    nrVerticesPerRectangle = 4
    nrLateralLinks = 4
    nrMaxLinks = 6

    latitude = NODATA
    longitude = NODATA
    elevation = NODATA          # [m] topographic elevation (center grid)
    timeZone = NODATA           # [-] time zone

    cellSize = NODATA           # [m] cell size
    gridWidth = NODATA          # [m] grid width (axis x)
    gridHeight = NODATA         # [m] grid height (axis y)
    gridDepth = NODATA          # [m] grid depth (axis z)

    slopeY = NODATA             # [-] slope on y axis
    plantSlope = NODATA         # [-] slope around plant (baulatura)
    plantSlopeWidth = NODATA    # [m] width of baulatura

    nrLayers = NODATA
    nrCells = NODATA
    totalArea = NODATA

    nrRectanglesInXAxis = NODATA
    nrRectanglesInYAxis = NODATA
    nrRectangles = NODATA


def initialize3DStructure(width, height, cellSize):
    C3DStructure.cellSize = cellSize

    C3DStructure.nrRectanglesInXAxis = int(width / cellSize)
    if (C3DStructure.nrRectanglesInXAxis % 2) == 0:
        C3DStructure.nrRectanglesInXAxis += 1

    C3DStructure.nrRectanglesInYAxis = int(height / cellSize)
    if (C3DStructure.nrRectanglesInYAxis % 2) == 0:
        C3DStructure.nrRectanglesInYAxis += 1

    C3DStructure.nrRectangles = C3DStructure.nrRectanglesInXAxis * C3DStructure.nrRectanglesInYAxis

    C3DStructure.gridWidth = C3DStructure.nrRectanglesInXAxis * cellSize
    C3DStructure.gridHeight = C3DStructure.nrRectanglesInYAxis * cellSize


class Clink:
    def __init__(self):
        self.index = NOLINK     # [-] index of linked cell
        self.area = NODATA      # [m^2] area of interface
        self.distance = NODATA  # [m]


class Cboundary:
    def __init__(self):
        self.type = BOUNDARY_NONE
        self.area = NODATA      # [m^2] area of interface
        self.slope = NODATA     # [-] slope
        self.flow = NODATA      # [m^3 s^-1] boundary water flow


class Ccell:
    def __init__(self):
        self.x = NODATA             # [m]
        self.y = NODATA             # [m]
        self.z = NODATA             # [m]
        self.volume = NODATA        # [m^3] volume
        self.area = NODATA          # [m^3] area (for surface cells)
        self.horizonIndex = NODATA  # soil horizon index
        self.isSurface = False      # true if the cell is on surface
        self.Se = NODATA            # [-] degree of saturation
        self.H = NODATA             # [m] current total water potential
        self.H0 = NODATA            # [m] total water potential
        self.k = NODATA             # [m s^-1] hydraulic conductivity
        self.sinkSource = NODATA    # [m^3 s^-1] water sink/source
        self.flow = NODATA          # [m^3 s^-1] sink/source + boundary
        self.boundary = Cboundary()
        self.upLink = Clink()
        self.downLink = Clink()
        self.lateralLink = [Clink(), Clink(), Clink(), Clink()]


# model parameters
class C3DParameters:
    # water retention curve and conductivity
    waterRetentionCurve = IPPISCH_VG
    conductivityMean = LOGARITHMIC
    conductivityHVRatio = 1.0

    # soil layers thickness
    minThickness = 0.01                 # [m]
    maxThickness = 0.04                 # [m]
    maxThicknessAt = 0.2                # [m]

    # processes
    computeInfiltration = True
    assignIrrigation = True
    assignPrecipitation = True
    computeEvaporation = True
    computeTranspiration = True

    # surface flow
    computeSurfaceFlow = False
    roughness = 0.24                    # [s m^0.33]
    pond = 0.005                        # [m]

    # boundary
    isSurfaceRunoff = False
    isFreeLateralDrainage = False
    isFreeDrainage = True
    isWaterTable = False

    # initial conditions
    initialWaterPotential = -3.0        # [m]
    waterTableDepth = -3.0              # [m]

    # numerical solution parameters
    currentDeltaT = 60.0                # [s]
    deltaT_min = 1                      # [s]
    deltaT_max = 3600.0                 # [s]
    currentDeltaT_max = deltaT_max      # [s]
    maxIterationsNr = 100
    maxApproximationsNr = 10
    residualTolerance = 1E-12
    MBRThreshold = 1E-5

    # simulation type
    isFirstAssimilation = True
    isPeriodicAssimilation = False
    isVisual = True
    assimilationInterval = 24


# global
C3DCells = []
dripperIndices = []
plantIndices = []
