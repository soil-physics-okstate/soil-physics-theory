#PSP_coupled1D.py
from __future__ import print_function, division
import numpy as np
import PSP_soil as soil
from PSP_public import *
import PSP_boundary as boundary  
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
import PSP_grid as grid

n = 20                                 

def initialize(mySoil, psi0, T0):
    global z, vol, dz, psi, oldPsi, dpsi, T, oldT
    global theta, oldTheta, vapour, oldVapour, Ch, oldCh, dHydr2
	
    z = grid.linear(n, mySoil.lowerDepth)
    dz = np.append(z[1:n+2]-z[0:n+1],z[n+1]-z[n])
    dz[0] = 1.
    vol = np.append(np.append(0,(z[2:n+1] - z[0:n-1]) / 2.0),(z[n] - z[n-2]) / 2.0)
    theta0 = soil.thetaFromPsi(soil.CAMPBELL, mySoil, psi0)
    Ch0 = soil.heatCapacity(mySoil, theta0) 
    psi = np.zeros(n+2,float)+psi0
    oldPsi = np.zeros(n+2,float)+psi0
    dPsi = np.zeros(n+2,float)
    T = np.zeros(n+2,float)+T0
    oldT = np.zeros(n+2,float)+T0
    theta = np.zeros(n+2,float)+theta0
    oldTheta = np.zeros(n+2,float)+theta0
    Ch = np.zeros(n+2,float)+Ch0
    oldCh = np.zeros(n+2,float)+Ch0
    dpsi = np.zeros(n+2,float)
	
def solver(mySoil, myBoundary, isFreeDrainage, dt):
    global psi, oldPsi, dpsi, T, oldT, theta
    global oldTheta, vapour, oldVapour, Ch, oldCh, dHydr2
	
    massBalance = 1.
    energyBalance = 1.
    nrIterations = 0
    psi = np.copy(oldPsi)
    theta = np.copy(oldTheta)
    T = np.copy(oldT)
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        T[n+1] = T[n]

    while (((massBalance > toleranceVapour) or (energyBalance > toleranceHeat)) 
											and (nrIterations < maxNrIterations)):
        Cw = soil.dTheta_dPsi(soil.CAMPBELL, mySoil, psi)
        meanPsi = np.append(0,(psi[2:n+2]+psi[1:n+1])/2.)  
        TCelsius = np.append(0,(T[2:n+2]+T[1:n+1])/2.) 
        gradPsi = np.append(0,(psi[2:n+2]-psi[1:n+1])/dz[1:n+1])
        gradT = np.append(0,(T[2:n+2]-T[1:n+1])/dz[1:n+1])
        FlowProp = soil.PropertiesForMassBalance(TCelsius[1:n+1],meanPsi[1:n+1],mySoil)
        dHydr = np.append(np.append(0,(theta[1:n+1]-oldTheta[1:n+1]) * vol[1:n+1]*1.e3 \
			+ (FlowProp.K11[0:n]*gradPsi[0:n]-FlowProp.K11[1:n+1]*gradPsi[1:n+1])*dt \
			+ (FlowProp.K21[0:n]*gradT[0:n]-FlowProp.K21[1:n+1]*gradT[1:n+1])*dt \
			+ (FlowProp.gravityWaterFlow[0:n]-FlowProp.gravityWaterFlow[1:n+1])*dt),0)
        aHydr = -np.append([0,0],FlowProp.K11[1:n+1]/dz[1:n+1]*dt)
        bHydr = np.append([0,0],FlowProp.K11[1:n+1]/dz[1:n+1]*dt)\
			+ np.append(np.append(0,Cw[1:n+1]*vol[1:n+1]*1e3),0)\
			+ np.append(np.append(0,FlowProp.K11[1:n+1]/dz[1:n+1]*dt),0)
        cHydr = -np.append(np.append(0,FlowProp.K11[1:n+1]/dz[1:n+1]*dt),0)
        airResistance = (1.0 / boundary.boundaryLayerConductance(myBoundary.windSpeed, myBoundary.airTemperature, T[1]))
        dHydr[1] -= boundary.waterFlux(psi[1],
        theta[1], T[1], myBoundary, airResistance) *dt
        bHydr[1] -= (boundary.dWaterFluxdPsi(psi[1], 
        theta[1], T[1], myBoundary, airResistance) *dt)
        massBalance = np.sum(abs(dHydr))
        ThomasBoundaryCondition(aHydr, bHydr,cHydr,dHydr, dpsi, 1, n)
        psi -= dpsi
        theta = soil.thetaFromPsi(soil.CAMPBELL, mySoil, psi)
        if (isFreeDrainage):
            psi[n+1] = psi[n]

        Ch = soil.heatCapacity(mySoil, theta)
        meanPsi = np.append(0,(psi[2:n+2]+psi[1:n+1])/2.)  
        TCelsius = np.append(0,(T[2:n+2]+T[1:n+1])/2.) 
        gradPsi = np.append(0,(psi[2:n+2]-psi[1:n+1])/dz[1:n+1])
        gradT = np.append(0,(T[2:n+2]-T[1:n+1])/dz[1:n+1])
        FlowProp = soil.PropertiesForMassBalance(TCelsius[1:n+1],meanPsi[1:n+1],mySoil)
        dTherm = np.append(np.append(0,Ch[1:n+1]*oldT[1:n+1]*vol[1:n+1] \
			+ (FlowProp.K12[0:n]*gradPsi[0:n]-FlowProp.K12[1:n+1]*gradPsi[1:n+1])*dt),0)
        aTherm = -np.append([0,0],FlowProp.K22[1:n+1]/dz[1:n+1]*dt)
        bTherm = np.append([0,0],FlowProp.K22[1:n+1]/dz[1:n+1]*dt) \
			+ np.append(np.append(0,Ch[1:n+1]*vol[1:n+1]),0) \
			+ np.append(np.append(0,FlowProp.K22[1:n]/dz[1:n]*dt),[0,0])
        cTherm = -np.append(np.append(0,FlowProp.K22[1:n]/dz[1:n]*dt),[0,0])
        airResistance = (1.0 / boundary.boundaryLayerConductance(myBoundary.windSpeed, myBoundary.airTemperature, T[1]))
        BoundaryThermalFlux = (boundary.thermalFlux(psi[1], theta[1], T[1], myBoundary, airResistance, False))
        dTherm[1] += BoundaryThermalFlux*dt
        bTherm[1] += (4.0*sigma*(myBoundary.airTemperature+273.15)**3*dt +(1200.0 / airResistance)*dt)
        ThomasBoundaryCondition(aTherm, bTherm, cTherm, dTherm, T, 1, n)
        
        #Energy balance for the surface layer
        energyBalance = (abs(Ch[1]*vol[1]*(T[1]-oldT[1]) +
        FlowProp.K22[1]/dz[1]*dt*(T[1]-T[2])\
			+ FlowProp.K12[1]*gradPsi[1]*dt \
			+ (4.0*sigma*(myBoundary.airTemperature+273.15)**3*dt
			+ (1200.0 / airResistance)*dt) * T[1] \
			- BoundaryThermalFlux*dt))
        #print(energyBalance)
        nrIterations += 1
		
    if ((massBalance < toleranceVapour) and (energyBalance < toleranceHeat)):
        isBalanceOk = True
        aerodynamicResistance = (1.0 / boundary.boundaryLayerConductance(
						myBoundary.windSpeed, myBoundary.airTemperature, T[1]))
        waterFlux = (boundary.waterFlux(psi[1], 
        theta[1], T[1], myBoundary, aerodynamicResistance))
        heatFlux = (boundary.thermalFlux(psi[1], 
        theta[1], T[1], myBoundary, aerodynamicResistance, True))
        evaporationFlux = boundary.evaporationFlux(psi[1], theta[1], T[1], myBoundary.airTemperature, 
						myBoundary.relativeHumidity, aerodynamicResistance)
        boundaryLayerConductance = boundary.boundaryLayerConductance(myBoundary.windSpeed, 
						myBoundary.airTemperature, T[1])
    else:
        isBalanceOk = False
        waterFlux = NODATA
        heatFlux = NODATA
        boundaryLayerConductance = NODATA
        evaporationFlux = NODATA
		
    return (isBalanceOk, waterFlux, heatFlux, boundaryLayerConductance, evaporationFlux, nrIterations, massBalance)
