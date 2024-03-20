import matplotlib.pyplot as plt
import PSP_thermalCond as soil
from PSP_readDataFile import *

def main():
    A, isFileOk = readDataFile("soilTemperatures_hours.csv", 1, ',', False)
    if (isFileOk == False): 
        print ("Incorrect format")
        return()
    
    soilTemperature = A[0]
    print ("Temperatures = ", soilTemperature)
    nrTemperatures = len(soilTemperature)
    
    myStr = "bulk density [kg/m^3]: " 
    bulkDensity = float(input(myStr))
    myStr = "clay [0 - 1]: " 
    clay = float(input(myStr))

    particleDensity = 2650
    porosity = 1 - (bulkDensity / particleDensity) 
    
    step = 0.02  
    nrValues = int(porosity / step) + 1
    waterContent = np.zeros(nrValues)
    thermalConductivity = np.zeros(nrValues)  
    for i in range(nrValues):
        waterContent[i] = step*i
    
   
    fig = plt.figure(figsize=(10,8))
    plt.xlabel('Water Content [m$^{3}$ m$^{-3}$]',fontsize=20,labelpad=8)
    plt.ylabel('Thermal Conductivity [W m$^{-1}$C$^{-1}$]',fontsize=20,labelpad=8)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=8)
    #plt.xlim(0, 0.7)
    for t in range(nrTemperatures):
        for i in range(nrValues):
            thermalConductivity[i] = soil.thermalConductivity(bulkDensity, 
                                    waterContent[i], clay, soilTemperature[t])
        if (t == 0): plt.plot(waterContent, thermalConductivity,'k')  
        if (t == 1): plt.plot(waterContent, thermalConductivity,'--k')   
        if (t == 2): plt.plot(waterContent, thermalConductivity,'-.k') 
        if (t == 3): plt.plot(waterContent, thermalConductivity,':k')
    plt.show() 
main()
