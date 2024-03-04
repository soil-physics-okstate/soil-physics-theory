#PSP_plot.py
import matplotlib.pyplot as plt
import PSP_public as public

def handle_close(evt):
    public.isQuitRequired = True

def plot_start(endTime):
    global myPlot
    plt.ion()
    fig, myPlot = plt.subplots(3,3, figsize=(16, 9.5), dpi=80)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
    myPlot[0,0].set_ylim(0, 20) 
    #myPlot[1,0].set_ylim(0,40)
    myPlot[2,0].set_ylim(0,100)
    myPlot[2,0].set_xlabel("Day of year",fontsize=14,labelpad=6)
    myPlot[0,0].set_ylabel("Precipitation [mm]",fontsize=14,labelpad=6)
    myPlot[1,0].set_ylabel("Air temperature [C]",fontsize=14,labelpad=6)
    myPlot[2,0].set_ylabel("Relative humidity [%]",fontsize=14,labelpad=6)
    myPlot[2,2].set_xlabel("Day of year",fontsize=14,labelpad=4)
    myPlot[0,2].set_ylabel("Bound. layer cond. [s m$^{-1}$]",fontsize=14,labelpad=6)
    myPlot[1,2].set_ylabel("Cum. evap. [kg m$^{-2}$]",fontsize=14,labelpad=6)
    myPlot[2,2].set_ylabel("Soil temperature [C]",fontsize=14,labelpad=6)
    fig.canvas.mpl_connect('close_event', handle_close)
    
def plot_variables(doy, z, theta, T, psi, boundaryLayerConductance, evaporation, prec, airT, rh):
    # Precipitation
    myPlot[0,0].plot(doy, prec,'k^')
    # Air temperature
    myPlot[1,0].plot(doy, airT,'k.')
    # Relative humidity
    myPlot[2,0].plot(doy, rh, 'k.')
    
    # Water content
    myPlot[0,1].clear()
    myPlot[0,1].set_xlim([0.,0.5])
    myPlot[0,1].set_xlabel("Water content [m$^3$ m$^{-3}$]",fontsize=14,labelpad=6)
    myPlot[0,1].set_ylabel("Depth [m]",fontsize=14,labelpad=6)
    myPlot[0,1].plot(theta[1:len(z)], -z[1:len(z)], 'ko')
    # Water potential
    myPlot[1,1].clear()
    myPlot[1,1].set_xlim([-2000,0])
    myPlot[1,1].set_xlabel("Water potential [J kg$^{-1}$]",fontsize=14,labelpad=6)
    myPlot[1,1].set_ylabel("Depth [m]",fontsize=14,labelpad=6)
    myPlot[1,1].plot(psi[1:len(z)], -z[1:len(z)], 'ko')
    # Soil Temperature
    myPlot[2,1].clear()
    myPlot[2,1].set_xlabel("Soil temperature [C]",fontsize=14,labelpad=6)
    myPlot[2,1].set_xlim([-10.,40.])
    myPlot[2,1].set_ylabel("Depth [m]",fontsize=14,labelpad=6)
    myPlot[2,1].plot(T[1:len(T)], -z[1:len(T)], 'ko')
    
    #TIME
    # Boundary Layer Conductance
    myPlot[0,2].plot(doy, boundaryLayerConductance,'k.')
    # Cumulative Evaporation 
    myPlot[1,2].plot(doy, evaporation,'k.')
    # Soil temperature
    myPlot[2,2].plot(doy, T[1],'k.')
    myPlot[2,2].plot(doy, T[2], 'k^')
    
    plt.pause(0.00001)

def plot_end():
    plt.ioff()
    plt.show()
