#PSP_boundaryLayerConductance
from __future__ import print_function, division

from PSP_readDataFile import readDataFile
from PSP_boundaryLayerConductance import *
from PSP_plot import *
    
def main():
    A, isFileOk = readDataFile('weather.dat', 1, '\t', False)
    if not isFileOk: 
        print ("Incorrect format in row: ", A)
        return()
    
    airT = A[:,1]
    windSpeed = A[:,4]
    
    atmPressure = 100
    TsoilK = 15 + zeroKelvin
    nrHours = 24*7            #one week
    
    plot_start(nrHours)
    for i in range(nrHours):
        K = Kh(windSpeed[i], airT[i], TsoilK, atmPressure)
        plot_variables(i, windSpeed[i], airT[i], K)
        
    plot_end() 
main()
