#PSP_soluteTransportAnalytical
import numpy as np
from scipy import special as sp
import matplotlib.pyplot as plt

def main(): 
    tmax = 11    
    nrHours = 10                                       
    n = 100                                                         
    D = 1.0 
    fw = 2.5        
    v = 5.0
    mass = 3.0
    z = np.zeros(n)
    conc = np.zeros(n)  
    sqrtPi =   np.sqrt(np.pi) 
    for i in range(n):
        z[i] = i

    plt.ion()
    conc[0] = 0.0
    for t in range(1, nrHours):
        for i in range(1, n):
            a = v / (sqrtPi*np.sqrt(D*t)) * np.exp(-((z[i] - v*t)**2) / (4.0*D*t))
            b = (v**2)/(2*D) * np.exp((v*z[i])/D) * sp.erfc((z[i]+v*t) / np.sqrt(4.0*D*t))
            conc[i]= mass/fw * (a - b)                  
        plt.clf()
        plt.xlabel('Concentration [g cm$^{-2}$]',fontsize=16,labelpad=8)
        plt.ylabel('Depth [cm]',fontsize=16,labelpad=8)
        plt.xlim(0, 1)
        plt.plot(conc, -z, 'k-')
        plt.pause(0.1)
         
    plt.ioff()
    plt.show
    plt.savefig("analytical.eps", transparent = True)
main()
