# boundaryConditions.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

from math import sqrt
from dataStructures import *


def updateBoundary(deltaT):
    for i in range(C3DStructure.nrCells):
        # Initialize
        C3DCells[i].flow = 0.0
        C3DCells[i].boundary.flow = 0.0

        # sink/source: precipitation, irrigation, evapotranspiration
        if (C3DCells[i].sinkSource != NODATA) and (C3DCells[i].sinkSource != 0.0):
            # [m3 s-1]
            C3DCells[i].flow = C3DCells[i].sinkSource

        if C3DCells[i].boundary.type != BOUNDARY_NONE:
            slope = C3DCells[i].boundary.slope
            meanH = (C3DCells[i].H + C3DCells[i].H0) * 0.5

            if C3DCells[i].boundary.type == BOUNDARY_RUNOFF:
                Hs = meanH - (C3DCells[i].z + C3DParameters.pond)
                if Hs > EPSILON_METER and slope > 0:
                    boundaryArea = C3DCells[i].boundary.area * Hs
                    maxFlow = (Hs * C3DCells[i].area) / deltaT
                    # Manning equation [m3 s-1]
                    flow = ((boundaryArea / C3DParameters.roughness) * (Hs ** (2./3.)) * sqrt(slope))
                    C3DCells[i].boundary.flow = -min(flow, maxFlow)

            elif C3DCells[i].boundary.type == BOUNDARY_FREELATERALDRAINAGE:
                k = C3DCells[i].k * C3DParameters.conductivityHVRatio
                C3DCells[i].boundary.flow = -k * C3DCells[i].boundary.area * slope

            elif C3DCells[i].boundary.type == BOUNDARY_FREEDRAINAGE:
                C3DCells[i].boundary.flow = -C3DCells[i].k * C3DCells[i].upLink.area

            elif C3DCells[i].boundary.type == BOUNDARY_PRESCRIBEDTOTALPOTENTIAL:
                dzy = C3DStructure.slopeY * C3DCells[i].y
                prescribedH = C3DStructure.elevation - dzy + C3DParameters.waterTableDepth
                dH = prescribedH - C3DCells[i].H
                C3DCells[i].boundary.flow = C3DCells[i].k * dH/0.1 * C3DCells[i].area

            C3DCells[i].flow += C3DCells[i].boundary.flow

            # check on water surface
            if C3DCells[i].isSurface and C3DCells[i].flow < 0:
                Hs = C3DCells[i].H - C3DCells[i].z
                maxFlow = (Hs * C3DCells[i].area) / deltaT
                if abs(C3DCells[i].flow) > maxFlow:
                    C3DCells[i].flow = -maxFlow
