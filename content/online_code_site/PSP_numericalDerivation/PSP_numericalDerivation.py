#PSP_numericalDerivation.py
from PSP_readDataFile import readDataFile
import matplotlib.pyplot as plt

def firstDerivative5Points(y):
    myDerivative = []
    for i in range(0,2):
        myDerivative.append(0)
    for i in range(2, len(y)-2):
        dY = (1./(12.)) * (y[i-2] - 8.*y[i-1] + 8.*y[i+1] - y[i+2])
        myDerivative.append(dY)
    for i in range(len(y)-2, len(y)):
        myDerivative.append(0)
    return(myDerivative)
            
def main():
    myHour = []
    myTemp = []
    A, isFileOk = readDataFile("airTemp.dat", 1, '\t', True)
    if (isFileOk == False): print ("Incorrect format in row: ", A)
    myHour = A[:,0]
    myTemp = A[:,1]
    myDerivative = firstDerivative5Points(myTemp)
    
    fig = plt.figure(figsize=(10,8))
    plt.plot(myHour, myTemp, 'ko')        
    plt.plot(myHour, myDerivative,'k')
    plt.xlim(xmin=0,xmax=25)
    plt.ylim(ymin=-5,ymax=40)
    plt.xlabel('Time [hour]',fontsize=20,labelpad=8)
    plt.ylabel('Temperature [C]',fontsize=20,labelpad=8)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
    plt.show()
    
main()