#PSP_columnWaterContent
from __future__ import print_function, division

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.image as image
import numpy as np
import PSP_soil as soil
from PSP_integration import qsimp

def main():
    choice = 0
    print (soil.CAMPBELL,' Campbell')
    print (soil.VAN_GENUCHTEN,' van Genuchten')
    while (choice < soil.CAMPBELL) or (choice > soil.VAN_GENUCHTEN):
        choice = float(input("Choose water retention curve: "))
        if (choice < soil.CAMPBELL) or (choice > soil.VAN_GENUCHTEN):
            print('wrong choice.')
    soil.waterRetentionCurve = choice
	
    waterTableDepth = 2.0                       
    nrValues = 100
    step = waterTableDepth/nrValues           
    
    integral = qsimp(soil.waterContent, -waterTableDepth*9.81, 0)
    totalWaterContent = integral / 9.81
    print ("\nTotal water content [m2/m2]:", totalWaterContent)
    
    x = np.zeros(nrValues+1, float)
    y = np.zeros(nrValues+1, float)
    for i in range(nrValues+1):
        y[i] = step*i                           
        psi = y[i] - waterTableDepth           
        psi *= 9.81                             
        x[i] = soil.waterContent(psi)
    
    img = image.imread('soilProfile.png')
    matplot = plt.imshow(img, extent=[0,1,waterTableDepth,0])   
    plt.plot(x, y, 'k', linewidth='3')
    plt.xlim(0, 1)
    plt.ylim(waterTableDepth, 0)
    plt.title('')
    plt.xlabel('water content [m$^2$ m$^{-2}$]')
    plt.ylabel('depth [m]')
    plt.show()
	
main()