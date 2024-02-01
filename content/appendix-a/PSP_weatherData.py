#PSP_weatherData.py
from PSP_read6VarFile import read6VarFile
import matplotlib.pyplot as plt

def  computeDaily(year, doy, hourMinute, prec, airT, rh, meanT, cumPrec):
    mySumT = airT[0]
    mySumP = prec[0]
    nValues = 1
    for i in range(1, len(airT)):
        if (int(hourMinute[i] / 100) != 0):
            mySumT += airT[i]
            mySumP += prec[i]
            nValues += 1
        else:
            meanT.append(mySumT / float(nValues))
            cumPrec.append(mySumP)
            nValues = 1
            mySumT = airT[i]
            mySumP += prec[i]
    #last day
    meanT.append(mySumT / float(nValues))
    cumPrec.append(mySumP)
      
def main():
    year=[]
    doy = []
    hourMinute = []
    prec=[]
    airT  = []
    rh=[]
    
    read6VarFile("weather.txt", 1, '\t', year, doy, hourMinute, prec, airT, rh, False)
    meanT = []
    cumPrec=[]
    computeDaily(year, doy, hourMinute, prec,airT, rh, meanT, cumPrec)

    myGraph = plt.figure()
    myAxisLeft = myGraph.add_subplot(111)
    myAxisRight = myAxisLeft.twinx()
    plt.title('')
    myAxisLeft.set_xlabel('Time [doy]',fontsize=16)
    myAxisLeft.set_ylabel('Temperature [C]',fontsize=16)
    myAxisRight.set_ylabel('Cumulative precipitation [mm]',fontsize=16)

    #=======================DAILY AVERAGE TEMP PLOT===================== 
    for h in range(0, len(meanT)):
        myAxisLeft.plot(h, meanT[h], 'ro')
        
    #=======================CUMULATIVE PRECIPITATION===================         
    h = range(0, len(cumPrec))
    myAxisRight.plot(h, cumPrec,'b-')

    plt.show() 
main()
