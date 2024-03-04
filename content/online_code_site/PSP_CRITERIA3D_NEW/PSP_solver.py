#PSP_solver.py
from PSP_waterProcesses import *
from PSP_public import C3DParameters, NOLINK, NODATA
from PSP_dataStructures import C3DCells, C3DStructure
import PSP_boundaryConditions as boundary
import numpy as np

A = np.array([[],[]], np.float64)
C = np.array([], np.float64)
x = np.array([], np.float64)
b = np.array([], np.float64)
indices = np.array([[],[]], int)

def setCriteria3DArrays(nrCells, nrLinks):
    global x, b, C, A, indices
    x.resize(nrCells)
    b.resize(nrCells)
    C.resize(nrCells)
    A.resize((nrCells, nrLinks))
    indices.resize((nrCells, nrLinks))
    
def computeStep(deltaT):
    global x, C, indices
    #initialize
    approximation = 1
    isValidStep = False
    for i in range(C3DStructure.nrCells):
        C3DCells[i].H0 = C3DCells[i].H
        x[i] = C3DCells[i].H
        if  (C3DCells[i].isSurface): 
            C[i] = C3DCells[i].area
    
    while ((not isValidStep) 
    and (approximation <= C3DParameters.maxApproximationsNr)):
        isFirstApprox = (approximation == 1)
        balance.maxCourant = 0.0
        for i in range(C3DStructure.nrCells):
            if (not C3DCells[i].isSurface):
                C3DCells[i].Se = soil.getDegreeOfSaturation(i)
                C3DCells[i].k = soil.getHydraulicConductivity(i)
                C[i] = C3DCells[i].volume * soil.getdTheta_dH(i)
        boundary.updateBoundary(deltaT)
        
        print ("approximation nr:", approximation)
        print ("Sum flows (abs) [m^3]:", format(balance.sumWaterFlow(deltaT, True),".3f"))  
        
        for i in range(C3DStructure.nrCells):
            k = 0
            if (newMatrixElement(i, C3DCells[i].upLink, k, 
                                 False, deltaT, isFirstApprox)): 
                k += 1
            for l in range(C3DStructure.nrLateralLinks):
                if (newMatrixElement(i, C3DCells[i].lateralLink[l], k,
                                 True, deltaT, isFirstApprox)): 
                    k += 1
            if (newMatrixElement(i, C3DCells[i].downLink, k,
                                 False, deltaT, isFirstApprox)): 
                k += 1
            if (k < C3DStructure.nrMaxLinks): 
                indices[i][k] = NOLINK
                
            arrangeMatrix(i, deltaT)
            
        if ((balance.maxCourant > 1.0) 
        and (deltaT > C3DParameters.deltaT_min)):
            print ("Courant too high:", balance.maxCourant)
            print ("Decrease time step")
            while (balance.maxCourant > 1.0):
                balance.halveTimeStep()
                balance.maxCourant *= 0.5
            return(False)

        if not solveMatrix(approximation):
            balance.halveTimeStep()
            print("System not convergent.")
            return(False)
        # check surface error
        for i in range(C3DStructure.nrTriangles):
            if (C3DCells[i].isSurface):
                if (x[i] < C3DCells[i].z):
                    x[i] = C3DCells[i].z
        # new hydraulic head
        for i in range(0, C3DStructure.nrCells):
            C3DCells[i].H = x[i]
            C3DCells[i].Se = soil.getDegreeOfSaturation(i)
        # balance
        isValidStep = balance.waterBalance(deltaT, approximation)
        if (balance.forceExit): return(False)
        approximation += 1
    
    return (isValidStep)

def  newMatrixElement(i, link, k, isLateral, deltaT, isFirstApprox):
    global A, indices
    j = link.index
    if (j == NOLINK): return (False)
    
    value = 0.0
    if C3DCells[i].isSurface:
        if C3DCells[j].isSurface:
            value = runoff(i, link, deltaT, isFirstApprox)
        else:
            value = infiltration(i, j, link, deltaT, isFirstApprox)
    else:
        if C3DCells[j].isSurface:
            value = infiltration(j, i, link, deltaT, isFirstApprox)
        else:
            value = redistribution(i, link, isLateral, deltaT)
    if (value == 0.0): return(False)
    
    indices[i][k] = j
    A[i][k] = value
    return (True)

def arrangeMatrix(i, deltaT):
    global b, C, A 
    mySum = 0.0
    for j in range(C3DStructure.nrMaxLinks):
        if (indices[i][j] == NOLINK): break             
        mySum += A[i][j]
        A[i][j] *= -1.0

    # diagonal and vector b     
    D = (C[i] / deltaT) + mySum
    b[i] = (C[i] / deltaT) * C3DCells[i].H0
    if (C3DCells[i].flow != NODATA):
        b[i] += C3DCells[i].flow
        
    # matrix conditioning
    b[i] /= D
    for j in range(C3DStructure.nrMaxLinks):
        if (indices[i][j] == NOLINK): break 
        A[i][j] /= D 

def solveMatrix(approximation):
    ratio = (C3DParameters.maxIterationsNr 
             / C3DParameters.maxApproximationsNr)
    maxIterationsNr = max(10, ratio * approximation)
    
    iteration = 0
    norm = 1000.
    bestNorm = norm
    while ((norm > C3DParameters.residualTolerance) 
    and (iteration < maxIterationsNr)): 
        norm = GaussSeidel()
        if norm > (bestNorm * 10.0): return(False)
        bestNorm = min(norm, bestNorm)
        iteration += 1
    return(True)

def GaussSeidel():
    global x
    norm = 0.0
    for i in range(C3DStructure.nrCells):
        new_x = b[i]
        for j in range(C3DStructure.nrMaxLinks):
            n = indices[i][j] 
            if (n == NOLINK): break
            new_x -= (A[i][j] * x[n])
                
        dx = fabs(new_x - x[i])
        # infinite norm
        if (dx > norm): norm = dx
        x[i] = new_x
    return(norm)
