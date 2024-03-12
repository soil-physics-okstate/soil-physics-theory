#PSP_plot.py
import matplotlib.pyplot as plt
    global myPlot
    f, myPlot = plt.subplots(3, figsize=(12,12))
    for i in range(3):
        myPlot[i].set_xlim(0, nrHours)
        
    myPlot[2].set_xlabel("Time [h]",fontsize=14,labelpad=6)
    myPlot[0].set_ylim(0, 6)
    myPlot[1].set_ylim(10, 30)
    myPlot[2].set_ylabel("Bound. layer cond. [s m$^{-1}$]",fontsize=14,labelpad=6)
    myPlot[2].set_ylim(0, 0.02)
def plot_variables(hour, windSpeed, airT, boundaryLayerConductance):
    plt.show()