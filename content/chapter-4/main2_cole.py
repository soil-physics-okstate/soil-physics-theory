import matplotlib.pyplot as plt
import numpy as np
from PSP_heat_mod import *  
from PSP_grid import *

def main(bulkDensity, clay):  
    global z
    print (FIN_DIFF, 'Finite Difference')
    print (CELL_CENT_FIN_VOL, 'Cell-Centered Finite Volume')
    
    solver = int(1)
    uniform = int(1)
    if (uniform == 1):
      myStr = "water content (m^3/m^3): " 
      theta = float(0.35)
      thetaIni = np.ones(22)*theta
    else:
      myStr = "surface water content (m^3/m^3): " 
      thetaSurf = float(input(myStr))
      myStr = "sub-surface water content (m^3/m^3): " 
      thetaSubSurf = float(input(myStr))
      thetaIni = np.concatenate((np.array([thetaSurf]), np.ones(21)*thetaSubSurf))
    # print(thetaIni)
    # print(len(thetaIni))

    myStr = "mean temperature [C]: " 
    meanT = float(25)
    myStr = "amplitude of change in temperature [C]: " 
    ampT = float(10)
    omega = 2.0 * np.pi/(24 * 3600.0)
    airT0 = meanT
    timeShift = 8                       

    if (solver == FIN_DIFF):
        myStr = "weighting factor for time discretization:"
        myStr += " (0: explicit, 1: implicit Euler) = " 
        factor = float(0.6)
 
    z = initialize(airT0, thetaIni, solver)
    simulationLength = int(48)    
                    
    endTime = simulationLength * 3600.0         
    timeStepMax = 3600.0                        
    dt = timeStepMax / 8.0                      
    time = 0.0                                  
    sumHeatFlux = 0
    totalIterationNr = 0
    
    f, plot = plt.subplots(3, figsize=(8,8), dpi=80)
    plt.subplots_adjust(hspace = 0.3)
    plot[0].set_xlabel("Temperature [C]",fontsize=14,labelpad=2)
    plot[0].set_ylabel("Depth [m]",fontsize=14,labelpad=4) 
    plot[0].set_xlim(meanT-ampT, meanT+ampT) 
    plot[1].set_xlabel("Time [h]",fontsize=14,labelpad=2)  
    plot[1].set_ylabel("Temperature [C]",fontsize=14,labelpad=4)
    plot[2].set_xlabel("Time [h]",fontsize=14,labelpad=2)
    plot[2].set_ylabel("Heat flux [W m$^{-2}$]",fontsize=14,labelpad=4)
    plot[1].set_xlim(timeShift, simulationLength+timeShift)
    plot[1].set_ylim(meanT-ampT, meanT+ampT)
    plot[2].set_xlim(timeShift, simulationLength+timeShift)
  
    #Create new output file in append mode
    import os
    try:
      os.remove('output.csv')
    except OSError:
      pass
    outFile= open("output.csv","a")

    #Write header with depth
    outFile.write("time [hr], 0.0 [m], 0.1 [m], 0.3 [m]\n")    

    while (time < endTime):
        dt = min(dt, endTime - time)
        airT = airT0 + ampT * np.sin((time+dt)*omega)
        if (solver == FIN_DIFF):
            success, nrIterations, heatFlux = (
                finiteDifference(airT, meanT, dt, factor, bulkDensity, clay))            
        elif (solver == CELL_CENT_FIN_VOL):
            success, nrIterations, heatFlux = (
                cellCentFiniteVol(airT, meanT, dt, bulkDensity, clay))
        totalIterationNr += nrIterations
        
        
        if success:
            #Convergence achieved
            for i in range(n+1):
                oldT[i] = T[i]
            sumHeatFlux += heatFlux * dt 
            time += dt
            
            t = time/3600. + timeShift
            
            plot[0].plot(T[1:len(T)], -z[1:len(T)], 'k')
            # plot[0].plot(T[1:len(T)], -z[1:len(T)], 'ko')
            plot[0].draw(plt.gcf().canvas.get_renderer())
            plot[1].plot(t, T[getLayerIndex(z, 0.0)], 'ko')    
            plot[1].plot(t, T[getLayerIndex(z, 0.1)], 'ks')     
            plot[1].plot(t, T[getLayerIndex(z, 0.3)], 'k^')
            plot[1].draw(plt.gcf().canvas.get_renderer())    
            plot[2].plot(t, heatFlux, 'ko')
            plot[2].draw(plt.gcf().canvas.get_renderer())
            #plt.pause(0.0001)
            #increment time step when system is converging
            if (float(nrIterations/maxNrIterations) < 0.25): 
                    dt = min(dt*2, timeStepMax)
            #print(heat.T[getLayerIndex(heat.z, 0.15)])
            #save to file, the number after z is the selected depth        
            outFile.write("%.3f,%.3f,%.3f,%.3f,\n" %(t,T[getLayerIndex(z, 0.0)],T[getLayerIndex(z, 0.1)],T[getLayerIndex(z, 0.3)]))

    
        else:
            #No convergence
            dt = max(dt / 2, 1)
            for i in range(n+1): T[i] = oldT[i]
            print ("dt =", dt, "No convergence")

    outFile.close()
    print("nr of iterations per hour:", totalIterationNr / simulationLength)
    #plt.ioff()
    plt.show()
# main()
