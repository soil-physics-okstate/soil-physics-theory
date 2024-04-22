from __future__ import print_function, division
from PSP_readDataFile import readDataFile
from PSP_TTwaterContent import *
from PSP_travelTime import *
import matplotlib.pyplot as plt


def importData(): 
    headerNr = 8 # number of header rows
    nrHeaderValues = int(headerNr)
    fileName = input("File name:")
    if (fileName != ""):
        global isDataLoaded, waveFormNrpoints
        x, isFileOk = readDataFile(fileName,nrHeaderValues,'\t', False)
        if (not isFileOk):
            showerror("Wrong file", "Wrong file!\nError reading row nr." + str(x))
            return(False)
        if len(x) == 1:
            data = x[0,:]
        else:
            data = x[:,0]

        reflecCoeff = normalizeVector(data)
        # print("number of values:", len(data))
        isDataLoaded = True
    return reflecCoeff

def drawRegressionLines():
    nrPoints = len(timeVector)
    step = int(16. * (nrPoints / 256.0))
    
    t = np.zeros(nrPoints, float)
    line1vec = np.zeros(nrPoints, float)
    line2vec = np.zeros(nrPoints, float)
    flatLinevec = np.zeros(nrPoints, float)
    for i in range(nrPoints): 
        t[i] = timeVector[i] * 1E09
        line1vec[i] = line1.a * timeVector[i] + line1.b
        line2vec[i] = line2.a * timeVector[i] + line2.b
        flatLinevec[i] = lastFlatLine.b
    
    index = int(p2.x / deltaTime)   
    first = max(0, index - step)
    last = min(nrPoints, index + step)
    plt.plot(t[first:last], line1vec[first:last], 'k')
    plt.plot(t[first:last], line2vec[first:last], 'k')
    
    first = nrPoints - step*4
    last = nrPoints
    plt.plot(t[first:last], flatLinevec[first:last], 'k')
    
    plt.plot(p0.x* 1E09, p0.y, 'rs')
    plt.plot(p1.x* 1E09, p1.y, 'rs')
    plt.plot(p2.x* 1E09, p2.y, 'rs')
    plt.plot(p3.x* 1E09, p3.y, 'rs')	
    
def main(waterTemperature,bulkDensity,solidPermittivity):
  
  isDataLoaded = False
  #set parameters
  #waterTemperature = int(input("Soil temperature (C):")) #20 C 
  #bulkDensity = int(input("Bulk density (kg/m^3):")) #1350 kg/m^3
  #solidPermittivity = int(input("Permittivity of soil solids (-): ")) #default is 4
  
  liquidPermittivity = getLiquidPermittivity(waterTemperature)
  headerNr = 8 # number of header rows
  vp = 0.99 	#fraction of speed of light [-]
  probeLength = 0.15 # m
  windowBegin = 0.
  windowWidth = 5.
  probeHandle = 0.12 # m
  handlePermittivity = 1.7
  point0X = 0
  point1X = 0
  point2X = 0
  v1 = 0
  vf = 0
  ratio = 0
  estBulkDensity = 0
  tt = 0
  bulkPermittivity = 0
  geometricPar = 0.5
  a = 0.0
  b = 0.14
  c = 1.2
  wcTopp = 0
  wcMalicki = 0
  wcMixModel = 0
      
  reflecCoeff = importData()
          
  if not isDataLoaded:
    showerror("Warning", "Data not loaded")
  
  #compute
  nrPoints = len(reflecCoeff)
  WF_parameters(vp, probeHandle, windowBegin, windowWidth, nrPoints) # calculates global variables deltaTime, deltaSpace, timeVector
  
  if not computeTravelTime(probeHandle, handlePermittivity, vp, reflecCoeff): #attempts to compute travel time and returns warning if error encountered
    showerror("Warning", "Wrong data, header or parameter")
  
  travelTime = p2.x - p1.x 
  bulkPermittivity = getBulkPermittivity(probeLength, travelTime, vp)
  wcTopp = getWaterContentTopp(bulkPermittivity)
  wcMalicki = getWaterContentMalicki(bulkPermittivity, bulkDensity)
  wcMixModel = getWaterContentMixModel(bulkPermittivity, bulkDensity, 
                      solidPermittivity, liquidPermittivity, geometricPar)
  
  print('Water content from mixing model')
  print(wcMixModel)
  return wcMixModel
  
  y1 = p1.y
  y2 = p2.y
  y3 = p3.y
  v1 = abs(y2-y1)
  vf = abs(y3+1)
  v1vf_ratio = v1/vf
  bd = int(getDensityCurioni(bulkPermittivity, v1, vf, a, b, c))
  
  #graph
  plt.close()
  plt.figure(figsize=(10,8))
  
  lastIndex = len(reflecCoeff)-2
  t = np.zeros(lastIndex, float)
  for i in range(lastIndex): 
    t[i] = timeVector[i] * 1E09
  y = reflecCoeff[0:lastIndex]
  dy = dy[0:lastIndex] 
  dy2 = dy2[0:lastIndex]   
  plt.plot(t, y, 'k.')
  #plt.plot(t, dy, 'k--')
  #plt.plot(t, dy2, 'r--')
  
  drawRegressionLines()
  
  plt.title("")
  plt.xlabel("Time [ns]",fontsize=20,labelpad=8)
  plt.ylabel("Reflection coefficient [-]",fontsize=20,labelpad=8)
  plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
  plt.tick_params(axis='both', which='minor', labelsize=20,pad=8)
  plt.ylim(bottom=-1, top=1)
  plt.show()
