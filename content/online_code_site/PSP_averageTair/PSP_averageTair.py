#PSP_averageTair.py
from PSP_read3VarFile import *
import matplotlib.pyplot as plt

def  computeMeanT(hourMinute, airT, meanT):
    mySum = 0. 
    myHour = 0
    nValues = 0
    for i in range(0, len(airT)):
        if (int(hourMinute[i] / 100) == myHour):
            mySum += airT[i]
            nValues += 1
        else:
            meanT.append(mySum / nValues)
            nValues = 1
            mySum = airT[i]
            myHour += 1
            if (myHour == 24): myHour = 0

def main():
    doy = []
    hourMinute = []
    airT  = []
    read3VarFile("tenMinutesTemp.txt", 1, '\t', doy, hourMinute, airT, False)
    meanT = []
    computeMeanT(hourMinute, airT, meanT) 
    
    for h in range(0, len(meanT)):
        plt.plot(h,meanT[h],'ro')
    plt.title('')
    plt.xlabel('Time [hour]',fontsize=16)
    plt.ylabel('Temperature [C]',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=14)

    plt.show()
main()
