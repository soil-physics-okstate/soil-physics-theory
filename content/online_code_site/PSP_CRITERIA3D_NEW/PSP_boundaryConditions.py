#PSP_boundaryConditions.py
from math import sqrt
from PSP_dataStructures import *
import PSP_soil as soil

def updateBoundary(deltaT):
    retentionCurve = C3DParameters.waterRetentionCurve 
    for i in range(C3DStructure.nrCells):
        #sink/source: precipitation, evapotranspiration
        if (C3DCells[i].sinkSource != NODATA):
            #[m3 s-1]
            C3DCells[i].flow = C3DCells[i].sinkSource
        else:      
            C3DCells[i].flow = 0.0
            
        if (C3DCells[i].boundary.type != BOUNDARY_NONE):
            C3DCells[i].boundary.flow = 0.0
            slope = C3DCells[i].boundary.slope 
            meanH = (C3DCells[i].H + C3DCells[i].H0) * 0.5;
                 
            if (C3DCells[i].boundary.type == BOUNDARY_RUNOFF):
                if (slope > 0.0):
                    Hs = meanH - (C3DCells[i].z + C3DParameters.pond)
                    if (Hs > EPSILON_METER):
                        boundaryArea = C3DCells[i].boundary.area * Hs
                        maxFlow = (Hs * C3DCells[i].area) / deltaT
                        # Manning equation [m3 s-1]
                        flow = ((boundaryArea / C3DParameters.roughness) 
                                            * (Hs**(2./3.)) * sqrt(slope))
                        C3DCells[i].boundary.flow = -min(flow, maxFlow)
                    
            elif (C3DCells[i].boundary.type == BOUNDARY_FREELATERALDRAINAGE): 
                if (slope > 0.0):              
                    signPsi = meanH - C3DCells[i].z
                    Se = soil.degreeOfSaturation(retentionCurve, signPsi)
                    k = soil.hydraulicConductivity(retentionCurve, Se)
                    k *= C3DParameters.conductivityHVRatio
                    C3DCells[i].boundary.flow = - k * C3DCells[i].boundary.area * slope
                                      
            elif (C3DCells[i].boundary.type == BOUNDARY_FREEDRAINAGE):
                signPsi = meanH - C3DCells[i].z
                Se = soil.degreeOfSaturation(retentionCurve, signPsi)
                k = soil.hydraulicConductivity(retentionCurve, Se)
                C3DCells[i].boundary.flow = -k * C3DCells[i].upLink.area
                
            C3DCells[i].flow += C3DCells[i].boundary.flow
            
            
            