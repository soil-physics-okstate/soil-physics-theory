#PSP_plot.py
import matplotlib.pyplot as pltdef plot_start(nrHours):
    global myPlot
    f, myPlot = plt.subplots(3, figsize=(12,12))    f.subplots_adjust(hspace=0.1)
    for i in range(3):
        myPlot[i].set_xlim(0, nrHours)
        
    myPlot[2].set_xlabel("Time [h]",fontsize=14,labelpad=6)    myPlot[0].set_ylabel("Wind speed [m/s]",fontsize=14,labelpad=6)
    myPlot[0].set_ylim(0, 6)    myPlot[1].set_ylabel("Air Temperature [C]",fontsize=14,labelpad=6)
    myPlot[1].set_ylim(10, 30)
    myPlot[2].set_ylabel("Bound. layer cond. [s m$^{-1}$]",fontsize=14,labelpad=6)
    myPlot[2].set_ylim(0, 0.02)    
def plot_variables(hour, windSpeed, airT, boundaryLayerConductance):    myPlot[0].plot(hour, windSpeed,'k.')    myPlot[1].plot(hour, airT,'ro')    myPlot[2].plot(hour, boundaryLayerConductance,'k*')def plot_end():
    plt.show()
