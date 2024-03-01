#PSP_waterProcesses.py
from math import fabs, sqrt
from PSP_public import C3DParameters, EPSILON_METER
from PSP_dataStructures import C3DCells
import PSP_soil as soil
import PSP_balance as balance

def redistribution(i, link, isLateral, deltaT):
    j = link.index
    k = soil.meanK(C3DParameters.meanType, C3DCells[i].k, C3DCells[j].k)
    if (isLateral):
        k *= C3DParameters.conductivityHVRatio
        
    return (k * link.area) / link.distance


def infiltration(surf, sub, link, deltaT, isFirstApprox):
    if (C3DCells[surf].z > C3DCells[sub].H):
        #unsaturated
        Havg = (C3DCells[surf].H + C3DCells[surf].H0) * 0.5
        Hs = Havg - C3DCells[surf].z
        if isFirstApprox:
            rain =  (C3DCells[surf].sinkSource 
                     / C3DCells[surf].area) * (deltaT / 2.0)
            Hs += rain
        if (Hs < 1E-12): return 0.0
        
        interfaceK = soil.meanK(C3DParameters.meanType, 
                                C3DCells[sub].k, soil.C3DSoil.Ks)
        dH = C3DCells[surf].H - C3DCells[sub].H
        maxK = (Hs / deltaT) * (link.distance / dH)
        k = min(interfaceK , maxK)
    else:
        #saturated
        k = soil.C3DSoil.Ks
    
    return (k  * link.area) / link.distance


def runoff(i, link, deltaT, isFirstApprox):
    j = link.index
    zmax = max(C3DCells[i].z, C3DCells[j].z)
    Hmax = max((C3DCells[i].H + C3DCells[i].H0)/ 2.0, 
               (C3DCells[j].H + C3DCells[j].H0)/ 2.0)
    Hs = Hmax - (zmax + C3DParameters.pond) 

    if isFirstApprox:
        rain = (C3DCells[i].sinkSource / C3DCells[i].area) * (deltaT / 2.0)
        Hs += rain
    if (Hs <= EPSILON_METER): return(0.0)
    
    dH = fabs(C3DCells[i].H - C3DCells[j].H)
    if (dH < EPSILON_METER): return (0.0)
    
    # [m/s] Manning equation
    v = (pow(Hs, 2.0 / 3.0) * sqrt(dH/link.distance)) / C3DParameters.roughness
    Courant = v * deltaT / link.distance
    balance.maxCourant = max(balance.maxCourant, Courant)

    # on surface: link.area = side length [m]
    area = link.area * Hs 
    return (v / dH) * area

