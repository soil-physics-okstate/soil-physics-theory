#PSP_coupled
from __future__ import print_function, division

import numpy as np
import PSP_public as public
import PSP_coupled1D as coupled
from PSP_soil import readSoil
from PSP_readDataFile import readGenericDataFile
from datetime import datetime
from PSP_plot import *
from PSP_longWaveRadiation import longWaveRadiationFromWeather
    
def main():
    isSuccess, mySoil = readSoil("soil.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
    A, isFileOk = readGenericDataFile('SPC_weather.txt', 1, '\t', False)
    if not isFileOk: 
        print ("Incorrect format in row: ", A)
        return()

    nrHours = len(A)
    strFirstDate = A[0][0]
    firstDate = datetime.strptime(strFirstDate, "%d/%m/%Y")
    firstDoy = firstDate.timetuple().tm_yday
    
    airT = np.zeros(nrHours, float)
    prec = np.zeros(nrHours, float)
    relativeHumidity = np.zeros(nrHours, float)
    windSpeed = np.zeros(nrHours, float)
    globalRadiation = np.zeros(nrHours, float)
    
    for i in range(nrHours):
        airT[i] = A[i][2]
        prec[i] = A[i][3]
        relativeHumidity[i] = A[i][4]
        windSpeed[i] = A[i][5]
        globalRadiation[i] = A[i][6]

    longWaveRadiation = longWaveRadiationFromWeather(firstDoy, airT, relativeHumidity, globalRadiation) 

    sumWaterFlow = 0.
    sumHeatFlow = 0.  
    sumEvaporationFlow = 0.  
    nrIterations = 0
    currentTime = 0                             #[s]
    dt = 300                                    #[s]
    oldDt = dt
    lastPlotTime = currentTime                  #[s]
    stepPlotTime = public.outputTimeStep*3600   #[s]
    endTime = 3600*nrHours                      #[s]
    
    coupled.initialize(mySoil, public.initialPotential, public.initialTemperature)
    if public.isPlotActivated: plot_start(endTime)
    
    outputFile = open(public.outputFileName, "w")
    outputFile.write('hour , water content [m^3 m^-3], water potential [J kg^-1], soil temperature [C]\n')
    layer = coupled.grid.getLayerIndex(public.outputDepth, coupled.z)
    
    while (currentTime < endTime) and (not public.isQuitRequired):
        h = int(currentTime/3600)
        startCurrentHour = h * 3600
        endCurrentHour = (h+1) * 3600
        if (currentTime == startCurrentHour): 
            dt = oldDt
        else: 
            oldDt = dt
            dt = min(dt, endCurrentHour - currentTime)

        myBoundary = coupled.boundary.Cboundary()
        myBoundary.time = currentTime
        myBoundary.airTemperature = airT[h]
        myBoundary.precipitation = prec[h]
        myBoundary.relativeHumidity = relativeHumidity[h]
        myBoundary.windSpeed = windSpeed[h]
        myBoundary.globalRadiation = globalRadiation[h]
        myBoundary.longWaveRadiation = longWaveRadiation[h]

        (isBalanceOk, waterFlux, heatFlux, boundaryLayerConductance, evaporationFlux, 
		nrIterations) = coupled.solver(mySoil, myBoundary, public.isFreeDrainage, dt)
        
        if isBalanceOk:
            for i in range(coupled.n+1):
                coupled.oldTheta[i] = coupled.theta[i]
                coupled.oldPsi[i] = coupled.psi[i]
                coupled.oldT[i] = coupled.T[i]
                coupled.oldCh[i] = coupled.Ch[i]
            sumWaterFlow += waterFlux * dt 
            sumHeatFlow += heatFlux * dt 
            sumEvaporationFlow += evaporationFlux * dt 
            currentTime += dt
            
            if  (currentTime - lastPlotTime) >= stepPlotTime:
                doy = firstDoy + currentTime/3600./24.
                if public.isPlotActivated: 
                    plot_variables(doy, coupled.z, coupled.theta, coupled.T, coupled.psi,
                        boundaryLayerConductance, sumEvaporationFlow, myBoundary.precipitation,
                        myBoundary.airTemperature, myBoundary.relativeHumidity)
                
                outputFile.write(str(h+1) 
                                 + ', ' + format(coupled.theta[layer],".2f") 
                                 + ', ' + format(coupled.psi[layer],".0f") 
                                 + ', ' + format(coupled.T[layer],".1f") 
                                 + "\n")
                lastPlotTime = currentTime
                
            if (float(nrIterations / public.maxNrIterations) < 0.1):  
                dt = int(min(dt*2, public.maxTimeStep))
                   
        else:
            print ("Hour =", h+1, "\tdt =", dt, "\tNo convergence")
            dt = max(int(dt / 2), 1)
            oldDt = dt
            for i in range(coupled.n+1):
                coupled.theta[i] = coupled.oldTheta[i]
                coupled.psi[i] = coupled.oldPsi[i]
                coupled.T[i] = coupled.oldT[i]
                
    if public.isPlotActivated and not public.isQuitRequired: plot_end()
    print("End!")
main()
