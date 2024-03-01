#PSP_coupled
from __future__ import print_function, division

import numpy as np
from PSP_public import *
import PSP_coupled1D as coupled
from PSP_soil import readSoil
from PSP_readDataFile import readGenericDataFile
from datetime import datetime
from PSP_plot import *
from PSP_longWaveRadiation import *
    
def main():
    isSuccess, mySoil = readSoil("soil.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
    A, isFileOk = readGenericDataFile('SPC_weather.txt', 1, '\t', False)
    if not isFileOk: 
        print ("Incorrect format in row: ", A)
        return()

    # read the data, convert from string to datetime and then extract
    # the DOY (day of year)
    nrHours = len(A)
    strFirstDate = A[0][0]
    firstDate = datetime.strptime(strFirstDate, "%d/%m/%Y")
    #print(firstDate)
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
    time = 0                        #[s]
    dt = 300                        #[s]
    lastPlotTime = time             #[s]
    stepHours = 1
    stepPlotTime = stepHours*3600   #[s]
    endTime = 3600*nrHours          #[s]
    
    coupled.initialize(mySoil, initialPotential, initialTemperature)
    plot_start(endTime)
    
    while (time < endTime):
        dt = min(dt, endTime - time)
        myBoundary = coupled.boundary.Cboundary()
        h = int(time/3600)
        myBoundary.time = time
        myBoundary.airTemperature = airT[h]
        myBoundary.precipitation = prec[h]
        myBoundary.relativeHumidity = relativeHumidity[h]
        myBoundary.windSpeed = windSpeed[h]
        myBoundary.globalRadiation = globalRadiation[h]
        myBoundary.longWaveRadiation = longWaveRadiation[h]

        (isBalanceOk, waterFlux, heatFlux, boundaryLayerConductance, evaporationFlux, 
		nrIterations,massBalance) = coupled.solver(mySoil, myBoundary, isFreeDrainage, dt)
 
        if isBalanceOk:
            for i in range(coupled.n+1):
                coupled.oldTheta[i] = coupled.theta[i]
                coupled.oldPsi[i] = coupled.psi[i]
                coupled.oldT[i] = coupled.T[i]
                coupled.oldCh[i] = coupled.Ch[i]
            sumWaterFlow += waterFlux * dt 
            sumHeatFlow += heatFlux * dt 
            sumEvaporationFlow += evaporationFlux * dt 
            time += dt
            
            if ((time - lastPlotTime) >= stepPlotTime):
                doy = firstDoy + time/3600./24.
                plot_variables(doy, coupled.z, coupled.theta, coupled.T, coupled.psi,
                        boundaryLayerConductance, sumEvaporationFlow, myBoundary.precipitation,
                        myBoundary.airTemperature, myBoundary.relativeHumidity)
                lastPlotTime = time
                
            if (float(nrIterations / maxNrIterations) <= 0.1):  
                dt = int(min(dt*2, maxTimeStep))
                    
        else:
            print ("time =", int(time), "\tdt =", dt, "\tNo convergence")
            dt = max(int(dt / 2), 1)
            for i in range(coupled.n+1):
                coupled.theta[i] = coupled.oldTheta[i]
                coupled.psi[i] = coupled.oldPsi[i]
                coupled.T[i] = coupled.oldT[i]
    plot_end()
main()
