import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('energy_balance.csv',delimiter=',',skiprows=1)
plt.plot(a[:,0],a[:,1],'ks', label='Longwave upward')
plt.plot(a[:,0],a[:,2],'-.k',label='Short wave downward')
plt.plot(a[:,0],a[:,3],'k^',label='Long wave downward')
plt.plot(a[:,0],a[:,4],'k',label='Latent heat')
plt.plot(a[:,0],a[:,5],'--k',label='Sensible heat')
plt.plot(a[:,0],a[:,6],':k',label='Soil heat flux')
plt.legend()
plt.xlabel('Time [h]',fontsize=20,labelpad=1)
plt.ylabel('Energy flux [W m$^{-2}$]',fontsize=20,labelpad=2)
plt.tick_params(axis='both', which='major', labelsize=16,pad=4)
plt.show()

