#PSP_readSoilTempData.py
from __future__ import print_function, division
from PSP_readDataFile import readGenericDataFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from datetime import datetime
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates

NODATA = -9999
            
def main():
    A, isFileOk = readGenericDataFile("soilTemperatures.txt", 1, ',', False)
    if (isFileOk == False): 
        print ("Incorrect format in row: ", A)
        return()
    
    nrValues = len(A)
    
    # date array
    date = []
    for i in range(nrValues):
        strDate = A[i][0]
        date.append(datetime.strptime(strDate,"%d-%m-%y %H.%M"))
    d = date2num(date)
    
    #values (with simple error check)
    soilTemp = np.zeros((4, nrValues))
    for col in range(4):
        previous = NODATA
        for i in range(nrValues):
            isOk = False
            # missing data
            if (A[i][col+1] != ""):
                value = float(A[i][col+1])
                # void data > 1 day
                if (i == 0) or ((d[i]-d[i-1]) < 1):
                    # wrong values
                    if (value > -50) and (value < 100):
                        # spikes
                        if (previous == NODATA) or (abs(value - previous) < 10):
                            isOk = True
            if isOk:
                soilTemp[col,i] = value
                previous = value
            else:    
                soilTemp[col,i] = np.nan
                previous = NODATA
    
    print ("first date:", date[0])
    print ("last date: ", date[len(date)-1])

    fig = plt.figure(figsize=(10,8))
    plt.title('')
    plt.ylabel('Temperature [C]', fontsize=20, labelpad=8)
    plt.xlabel('Time [date]', fontsize=20, labelpad=8)
    plt.plot_date(d, soilTemp[0],'k--')
    #plt.plot_date(d, soilTemp[1],'k--')
    plt.plot_date(d, soilTemp[2],'k-')
    plt.plot_date(d, soilTemp[3],'k.')
    plt.xticks(rotation='0')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%m.%y"))
    plt.show()
    
main()
