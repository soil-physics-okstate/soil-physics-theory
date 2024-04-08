from __future__ import print_function, division

from PSP_readDataFile import readDataFile
from PSP_TTwaterContent import *
from PSP_TTplot import *
from PSP_travelTime import computeTravelTime

import sys
if sys.version_info >= (3,0):
    from tkinter import *                           #3.x
    from tkinter.filedialog import askopenfilename      
    from tkinter.messagebox import showerror
else:
    from Tkinter import *                           #2.7
    from tkFileDialog import askopenfilename           
    from tkMessageBox import showerror
    
isDataLoaded = False
waterTemperature = 25
liquidPermittivity = getLiquidPermittivity(waterTemperature)

mainWindow = Tk()
mainWindow.title("TDRPy")
mainWindow.geometry("%dx%d" % (520, 700))

headerNrStr = StringVar()
vpStr = StringVar()             
probeLenghtStr = StringVar()
windowBeginStr = StringVar()
windowWidthStr = StringVar()
probleHandleStr = StringVar()
handlePermittivityStr = StringVar()
point0XStr = StringVar()
point1XStr = StringVar()
point2XStr = StringVar()
v1Str = StringVar()
vfStr = StringVar()
ratioStr = StringVar()
estBulkDensityStr = StringVar()
ttStr = StringVar()
bulkPermittivityStr = StringVar()
waterTemperatureStr = StringVar()
bulkDensityStr = StringVar()
aStr = StringVar()
bStr = StringVar()
cStr = StringVar()
solidPermittivityStr = StringVar()
liquidPermittivityStr = StringVar()
geometricParStr = StringVar()
wcToppStr = StringVar()
wcMalickiStr = StringVar()
wcMixModelStr = StringVar()


if sys.version_info >= (3, 0):
    waterTempLabel = Label(mainWindow, text="Water temp. [\u00B0C]")
else:
    waterTempLabel = Label(mainWindow, text="Water temp. [C]")
waterTempLabel.place(x=80, y=275)
waterTempWidget = Entry(mainWindow, width = 6, textvariable = waterTemperatureStr)
waterTempWidget.place(x=200, y=275)

def getEpsilonLabel():
    if sys.version_info >= (3, 0):
        return Label(mainWindow, font = "helvetica 12", text="\u03F5")
    else:
        return Label(mainWindow, font = "helvetica 10", text="e")
    
def checkTWater(event):
    global waterTemperature, liquidPermittivity
    if (event.widget == waterTempWidget):
        waterTemperature = float(waterTemperatureStr.get())
        liquidPermittivity = getLiquidPermittivity(waterTemperature)
        liquidPermittivityStr.set(format(liquidPermittivity,".3f"))

def importData(): 
    nrHeaderValues = int(headerNrStr.get())
    fileName = askopenfilename()
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
  breakpoint()
        tt.reflecCoeff = tt.normalizeVector(data)
        print("number of values:", len(data))
        isDataLoaded = True
        
def ComputeTT():
    if not isDataLoaded:
        showerror("Warning", "Data not loaded")
        return
    #read parameters
    vp = float(vpStr.get())                     
	#fraction of speed of light [-]
    probleLenght = float(probeLenghtStr.get())
    windowBegin = float(windowBeginStr.get())
    windowWidth = float(windowWidthStr.get())
    probeHandle = float(probleHandleStr.get())
    handlePermittivity = float(handlePermittivityStr.get())
    bulkDensity = float(bulkDensityStr.get())
    solidPermittivity = float(solidPermittivityStr.get())
    geomParameter = float(geometricParStr.get())
    a = float(aStr.get())
    b = float(bStr.get())
    c = float(cStr.get())
    
    #compute
    nrPoints = len(tt.reflecCoeff)
    tt.WF_parameters(vp, probeHandle, windowBegin, windowWidth, nrPoints)

    if not computeTravelTime(probeHandle, handlePermittivity, vp):
        showerror("Warning", "Wrong data, header or parameter")
        return
    
    travelTime = tt.p2.x - tt.p1.x
    bulkPermittivity = getBulkPermittivity(probleLenght, travelTime, vp)
    wcTopp = getWaterContentTopp(bulkPermittivity)
    wcMalicki = getWaterContentMalicki(bulkPermittivity, bulkDensity)
    wcMixModel = getWaterContentMixModel(bulkPermittivity, bulkDensity, 
                        solidPermittivity, liquidPermittivity, geomParameter)
    y1 = tt.p1.y
    y2 = tt.p2.y
    y3 = tt.p3.y
    v1 = abs(y2-y1)
    vf = abs(y3+1)
    v1vf_ratio = v1/vf
    bd = int(getDensityCurioni(bulkPermittivity, v1, vf, a, b, c))
    
    #print results
    x0 = tt.p0.x * (1E09)
    point0XStr.set(format(x0, '.3f'))
    x1 = tt.p1.x * (1E09)
    point1XStr.set(format(x1, '.3f'))
    x2 = tt.p2.x * (1E09)
    point2XStr.set(format(x2, '.3f'))
    ttStr.set(format(travelTime * (1E09), '.3f'))
    bulkPermittivityStr.set(format(bulkPermittivity, '.2f'))
    
    v1Str.set(format(v1, '.3f'))
    vfStr.set(format(vf, '.3f'))
    ratioStr.set(format(v1vf_ratio, '.3f'))
    estBulkDensityStr.set(bd)
    
    wcToppStr.set(format(wcTopp,".3f"))
    wcMalickiStr.set(format(wcMalicki,".3f"))
    wcMixModelStr.set(format(wcMixModel,".3f"))
    
    #graph
    cleanDisplay()
    drawWaveForm()
    drawRegressionLines()
    showDisplay()
      
def main():        
    vpStr.set(0.99)
    probeLenghtStr.set(0.15)
    windowBeginStr.set(0.)
    windowWidthStr.set(5.)
    probleHandleStr.set(0.120)
    handlePermittivityStr.set(1.7)
    point0XStr.set(0)
    point1XStr.set(0)
    point2XStr.set(0)
    
    v1Str.set(0)
    vfStr.set(0)
    ratioStr.set(0)
    estBulkDensityStr.set(0)
    ttStr.set(0)
    bulkPermittivityStr.set(0)
    bulkDensityStr.set(1350)
    waterTemperatureStr.set(waterTemperature)
    solidPermittivityStr.set(4.0)
    liquidPermittivityStr.set(format(liquidPermittivity, ".3f"))
    geometricParStr.set(0.5)
    aStr.set(0.0)
    bStr.set(0.14)
    cStr.set(1.2)
    wcToppStr.set(0)
    wcMalickiStr.set(0)
    wcMixModelStr.set(0)
    
    buttonImport = Button(mainWindow, text="Import data", command=importData)
    buttonImport.place(x=5, y=15)
    
    headerLabel = Label(mainWindow, text="Header values nr.")
    headerLabel.place(x=100, y=20)
    headerWidget = Entry(mainWindow, width = 3, textvariable = headerNrStr)
    headerWidget.place(x=200, y=20)
    headerWidget.insert(0, "8")
    headerWidget.pack

    computeTTButton = Button(mainWindow, text="Compute", command=ComputeTT)
    computeTTButton.place(x=250, y=15)
            
    vpLabel = Label(mainWindow, text="Vp [-]")
    vpLabel.place(x=80, y=80)
    vpWidget = Entry(mainWindow, width = 6, textvariable = vpStr)
    vpWidget.place(x=200, y=80)
    
    probeLenghtLabel = Label(mainWindow, text="Probe length [m]")
    probeLenghtLabel.place(x=80, y=105)
    probeLenghtWidget = Entry(mainWindow, width = 6, textvariable = probeLenghtStr)
    probeLenghtWidget.place(x=200, y=105)
    
    winBeginLabel = Label(mainWindow, text="Window begin [m]")
    winBeginLabel.place(x=80, y=130)
    winBeginWidget = Entry(mainWindow, width = 6, textvariable = windowBeginStr)
    winBeginWidget.place(x=200, y=130)
    
    winwidthLabel = Label(mainWindow, text="Window width [m]")
    winwidthLabel.place(x=80, y=155)
    winwidthWidget = Entry(mainWindow, width = 6, textvariable = windowWidthStr)
    winwidthWidget.place(x=200, y=155)
    
    probeHandleLabel = Label(mainWindow, text="Probe handle [m]")
    probeHandleLabel.place(x=80, y=180)
    probeHandleWidget = Entry(mainWindow, width = 6, textvariable = probleHandleStr)
    probeHandleWidget.place(x=200, y=180)
    
    epsilonLabel = getEpsilonLabel()
    epsilonLabel.place(x=80, y=200)
    epsilonLabel = Label(mainWindow, text="handle")
    epsilonLabel.place(x=90, y=206)
    permittivityWidget = Entry(mainWindow, width = 6, textvariable = handlePermittivityStr)
    permittivityWidget.place(x=200, y=200)
    
    # Soil parameters
    SoilParameterLabel = Label(mainWindow, text="Soil parameters" , font = "helvetica 10 bold", fg="red")
    SoilParameterLabel.place(x=20, y=225)
    
    bulkDensityLabel = Label(mainWindow, text="Bulk density [kg m-3]")
    bulkDensityLabel.place(x=80, y=250)
    bulkDensityWidget = Entry(mainWindow, width = 6, textvariable = bulkDensityStr)
    bulkDensityWidget.place(x=200, y=250)
    
    epsilon2Label = getEpsilonLabel()
    epsilon2Label.place(x=80, y=295)
    epsilon2Label = Label(mainWindow, text="liquid")
    epsilon2Label.place(x=90, y=300)
    liquidPermittivityLabel = Label(mainWindow, textvariable = liquidPermittivityStr)
    liquidPermittivityLabel.place(x=200, y=300)
    
    epsilon3Label = getEpsilonLabel()
    epsilon3Label.place(x=80, y=320)
    epsilon3Label = Label(mainWindow, text="solid")
    epsilon3Label.place(x=90, y=325)
    solidPermittivityWidget = Entry(mainWindow, width = 6, textvariable = solidPermittivityStr)
    solidPermittivityWidget.place(x=200, y=325)
    
    geometricParLabel = Label(mainWindow, text="alpha (geom. param.)")
    geometricParLabel.place(x=80, y=350)
    geometricParWidget = Entry(mainWindow, width = 6, textvariable = geometricParStr)
    geometricParWidget.place(x=200, y=350)
    
    # Bulk density parameters
    posY_col2 = 225
    posX_col2 = 320
    bdParamLabel = Label(mainWindow, text="Density est. parameters" , font = "helvetica 10 bold", fg="red")
    bdParamLabel.place(x=300, y=posY_col2)
    posY_col2 += 25
    
    aLabel = Label(mainWindow, text="a [-]")
    aLabel.place(x=posX_col2, y=posY_col2)
    aWidget = Entry(mainWindow, width = 6, textvariable = aStr)
    aWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25
    
    bLabel = Label(mainWindow, text="b [-]")
    bLabel.place(x=posX_col2, y=posY_col2)
    bWidget = Entry(mainWindow, width = 6, textvariable = bStr)
    bWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25
    
    cLabel = Label(mainWindow, text="c [-]")
    cLabel.place(x=posX_col2, y=posY_col2)
    cWidget = Entry(mainWindow, width = 6, textvariable = cStr)
    cWidget.place(x=posX_col2 + 120, y=posY_col2)
    posY_col2 += 25
    
    posY_col1 = 375
    posY_col2 = 375
    posX_col2 = 320
    
    # travel time output 
    ttResultsLabel = Label(mainWindow, text="Travel Time results" , font = "helvetica 10 bold", fg="blue")
    ttResultsLabel.place(x=20, y=posY_col1)
    posY_col1 += 25
    
    point0Label = Label(mainWindow, text="point 0 x [ns]")
    point0Label.place(x=80, y=posY_col1)
    point0Widget = Label(mainWindow, width = 6, textvariable = point0XStr)
    point0Widget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    point1Label = Label(mainWindow, text="point 1 x [ns]")
    point1Label.place(x=80, y=posY_col1)
    point1Widget = Label(mainWindow, width = 6, textvariable = point1XStr)
    point1Widget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    point2Label = Label(mainWindow, text="point 2 x [ns]")
    point2Label.place(x=80, y=posY_col1)
    point2Widget = Label(mainWindow, width = 6, textvariable = point2XStr)
    point2Widget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    ttLabel = Label(mainWindow, text="Travel Time [ns]")
    ttLabel.place(x=80, y=posY_col1)
    ttWidget = Label(mainWindow, width = 6, textvariable = ttStr)
    ttWidget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    bulkPermittivityLabel = Label(mainWindow, text="Bulk permittivity")
    bulkPermittivityLabel.place(x=80, y=posY_col1)
    bulkPermittivityWidget = Label(mainWindow, width = 6, textvariable = bulkPermittivityStr)
    bulkPermittivityWidget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    # Water content output
    wcLabel = Label(mainWindow, text="Water Content" , font = "helvetica 10 bold", fg="blue")
    wcLabel.place(x=20, y=posY_col1)
    posY_col1 += 25
    
    ToppLabel = Label(mainWindow, text="Topp")
    ToppLabel.place(x=80, y=posY_col1)
    ToppWidget = Label(mainWindow, width = 6, textvariable = wcToppStr)
    ToppWidget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    MalickiLabel = Label(mainWindow, text="Malicki")
    MalickiLabel.place(x=80, y=posY_col1)
    MalickiWidget = Label(mainWindow, width = 6, textvariable = wcMalickiStr)
    MalickiWidget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    dielecMixModelLabel = Label(mainWindow, text="Diel. mix model")
    dielecMixModelLabel.place(x=80, y=posY_col1)
    dielecMixModelWidget = Label(mainWindow, width = 6, textvariable = wcMixModelStr)
    dielecMixModelWidget.place(x=200, y=posY_col1)
    posY_col1 += 25
    
    # Density output
    reflectionLabel = Label(mainWindow, text="Density" , font = "helvetica 10 bold", fg="blue")
    reflectionLabel.place(x=300, y=posY_col2)
    posY_col2 += 25
    
    v1Label = Label(mainWindow, text="V1 [-]")
    v1Label.place(x=posX_col2, y=posY_col2)
    v1Widget = Label(mainWindow, width = 6, textvariable = v1Str)
    v1Widget.place(x=posX_col2 + 140, y=posY_col2)
    posY_col2 += 25
    
    vfLabel = Label(mainWindow, text="Vf [-]")
    vfLabel.place(x=posX_col2, y=posY_col2)
    vfWidget = Label(mainWindow, width = 6, textvariable = vfStr)
    vfWidget.place(x=posX_col2 + 140, y=posY_col2)
    posY_col2 += 25
    
    ratioLabel = Label(mainWindow, text="V1 / Vf [-]")
    ratioLabel.place(x=posX_col2, y=posY_col2)
    ratioWidget = Label(mainWindow, width = 6, textvariable = ratioStr)
    ratioWidget.place(x=posX_col2 + 140, y=posY_col2)
    posY_col2 += 25
    
    estBulkDensityLabel = Label(mainWindow, text="Density [m3 kg]")
    estBulkDensityLabel.place(x=posX_col2, y=posY_col2)
    estBulkDensityWidget = Label(mainWindow, width = 6, textvariable = estBulkDensityStr)
    estBulkDensityWidget.place(x=posX_col2 + 140, y=posY_col2)
    posY_col2 += 25

    mainWindow.bind("<Leave>", checkTWater)
    mainWindow.protocol("WM_DELETE_WINDOW", mainWindow.destroy)
        
    mainWindow.mainloop()
main()
    
