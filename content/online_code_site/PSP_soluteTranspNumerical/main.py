#PSP_main.py
from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import PSP_infiltration1D as inf
import PSP_diffusion1D as dif
    
def main():  
    isSuccess, soil = inf.readSoil("soil.txt")
    if not isSuccess: 
        print("warning: wrong soil file.")
        return
    
    funcType = inf.CAMPBELL
    solver = inf.NEWTON_RAPHSON_MFP
    thetaIni = 0.2
    #ubPotential = inf.airEntryPotential(funcType, soil) 
    ubPotential = inf.waterPotential(funcType,soil,thetaIni)  
   

    inf.initializeWater(funcType, soil, thetaIni, solver)
    dif.initializeDiffusion(inf.n, inf.z, inf.dz, inf.vol, inf.theta)
        
    
    print()
    print ("1: Free drainage")
    print ("2: Constant water potential")
    boundary = int(input("Select lower boundary condition:"))
    if (boundary == 1):
        isFreeDrainage = True
    else:
        isFreeDrainage = False
    
    simulationLenght = int(input("\nNr of simulation hours:"))    
                    
    endTime = simulationLenght * 3600   
    maxTimeStep = 3600                  
    dt = maxTimeStep / 10               
    time = 0                            
    sumInfiltration = 0
    totalIterationNr = 0
    
        
    plt.ion()
    f, myPlot = plt.subplots(3, figsize=(10, 8))
    f.subplots_adjust(hspace=.4)
    myPlot[0].set_xlim(0, 0.5)     
    myPlot[0].tick_params(axis='both', which='major', labelsize=14,pad=6)
    myPlot[1].set_xlim(0, simulationLenght * 3600)
    myPlot[1].set_ylabel("Infiltration Rate [kg m$^{-2}$ s$^{-1}$]",fontsize=14,labelpad=6)
    myPlot[1].set_xlabel("Time [s]",fontsize=14,labelpad=6)
    myPlot[1].tick_params(axis='both', which='major', labelsize=14,pad=6)
    myPlot[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    while (time < endTime):
        dt = min(dt, endTime - time)
        success, nrIterations, flux = inf.NewtonRapsonMFP(funcType,
                        soil, dt, ubPotential, isFreeDrainage)
        totalIterationNr += nrIterations
        
        if success:
            successDif = dif.cellCentFiniteVolWater(dt, inf.waterFluxDensity, inf.theta, inf.LOGARITHMIC, 1.)
            if not(successDif):
                 print("\nDiffusions solver does not converge. \n")
                 return
            for i in range(inf.n+2):
                inf.oldTheta[i] = inf.theta[i]
                dif.oldConc[i] = dif.concWat[i] * ( dif.klQ + (inf.waterDensity*inf.theta[i])/dif.bulkSoilDensity)
            sumInfiltration += flux * dt 
            time += dt
                        
            print("time =", int(time), "\tdt =", int(dt), 
                  "\tIter. Water =", int(nrIterations), 
                  "\tInf:", format(sumInfiltration, '.3f'))
            myPlot[0].clear()
            myPlot[0].set_xlabel("Water content [m$^3$ m$^{-3}$]",fontsize=14,labelpad=6)
            myPlot[0].set_ylabel("Depth [m]",fontsize=14,labelpad=6)
            myPlot[0].set_xlim([0,0.4])
            myPlot[0].plot(inf.theta[1:np.size(inf.theta)], -inf.z[1:np.size(inf.z)],'r-')
            myPlot[0].plot(inf.theta[1:np.size(inf.theta)], -inf.z[1:np.size(inf.z)], 'ko')
            myPlot[1].plot(time, flux, 'ko')
            myPlot[1].set_xlabel("Time [s]",fontsize=14,labelpad=6)
            myPlot[2].clear()
            myPlot[2].set_xlabel("Mass of solute per mass of soil [kg kg$^{-1}$]",fontsize=14,labelpad=6)
            myPlot[2].set_ylabel("Depth [m]",fontsize=14,labelpad=6)
            myPlot[2].plot(dif.oldConc[1:inf.n], -inf.z[1:inf.n], 'k')
            myPlot[2].set_xlim([0,2e-4])
            myPlot[2].tick_params(axis='both', which='major', labelsize=16,pad=6)
            myPlot[2].legend()
            plt.draw()
            
            if (float(nrIterations/inf.maxNrIterations) < 0.25): 
                    dt = min(dt*2, maxTimeStep)
        
        else:
            dt = max(dt / 2, 1)
            
            for i in range(inf.n+1):
                inf.theta[i] = inf.oldTheta[i]
                inf.psi[i] = inf.MFPFromTheta(soil, inf.theta[i])
            print ("dt =", dt, "No convergence")
            
    print("nr of iterations per hour:", totalIterationNr / simulationLenght)
    plt.ioff()
    plt.show()
main()
