#PSP_visual3D.py
from PSP_public import EPSILON_METER, C3DParameters
from PSP_dataStructures import C3DCells, C3DStructure
from PSP_color import *
from PSP_fileUtilities import loadState, saveState
import PSP_balance as balance
import PSP_tin as tin
import PSP_soil as soil
import vpython as visual

surfaceTriangles = []
subSurfaceTriangles = []
visualizedLayer = 1
nrColorLevels = 10
degreeMaximum = 1
degreeMinimum = 0.5
isPause = False
  

def initialize(totalWidth):
    global surfaceCanvas, soilCanvas, interface
    global layerLabel, timeLabel, precLabel, stepLabel, storageLabel
    global flowLabel, totalFlowLabel, totalErrorLabel, colorScale
    global visualizedLayer
    
    setAllColorScale()
    
    #CANVAS DIMENSION
    interfaceWidth = int(totalWidth * 0.2)
    dx = int((totalWidth-interfaceWidth) / 2.0)
    dy = int(dx * 0.8)
    h = int(dy / 30)
    
    #CENTER
    cX = (tin.header.xMin + tin.header.xMax) * 0.5
    cY = (tin.header.yMin + tin.header.yMax) * 0.5
    cZ = tin.header.zMin * tin.header.magnify
    Zlabel = (tin.header.zMax + tin.header.dz*0.75) * tin.header.magnify
    
    #INTERFACE CANVAS
    interface = visual.canvas(width = interfaceWidth, height = dy, align="left")
    interface.background = visual.color.white
    if (not C3DParameters.computeOnlySurface):
        interface.center = visual.vector(1,0,0)
        interface.range = 5
    else:
        interface.center = visual.vector(0,0,0)
        interface.range = 4

    timeLabel = visual.label(height=h, pos=visual.vector(0,3,0), text = "", canvas = interface)
    stepLabel = visual.label(height=h, pos=visual.vector(0,2,0), text = "", canvas = interface)
    precLabel = visual.label(height=h,pos=visual.vector(0,1,0), text = "", canvas = interface)
    storageLabel = visual.label(height=h,pos=visual.vector(0,0,0), text = "", canvas = interface)
    flowLabel = visual.label(height=h, pos=visual.vector(0,-1,0), text = "", canvas = interface)
    totalFlowLabel = visual.label(height=h, pos=visual.vector(0,-2,0), text = "", canvas = interface)
    totalErrorLabel = visual.label(height=h, pos=visual.vector(0,-3,0), text = "", canvas = interface)
    
    if (not C3DParameters.computeOnlySurface):
        #COLOR LEGEND
        stepY = 10 / nrColorLevels
        colorScale = []
        for i in range (nrColorLevels+1):
            label = visual.label(canvas = interface, pos=visual.vector(4, -5+(i*stepY), 0), height=h, 
                                 background = visual.vec(0,0,0))
            colorScale.append(label)

        #SOIL CANVAS
        soilCanvas = visual.canvas(width = dx, height = dy, align="left")
        soilCanvas.background = visual.color.white
        soilCanvas.center = visual.vector(cX, cY, cZ)
        soilCanvas.ambient = visual.vector(0.33, 0.33, 0.5)
        soilCanvas.up = visual.vector(0,0,1)
        soilCanvas.forward = visual.vector(0.33, -0.33, -0.15)
        soilCanvas.range = (tin.header.xMax - tin.header.xMin) * 0.55
        layerLabel = visual.label(canvas = soilCanvas, height = h, pos=visual.vector(cX, cY, Zlabel))
        
        drawColorScale()
        drawSubSurface(True)
    
    #SURFACE CANVAS
    surfaceCanvas = visual.canvas(width = dx, height = dy, align="left")
    surfaceCanvas.background = visual.color.white
    surfaceCanvas.center = visual.vector(cX, cY, cZ)
    surfaceCanvas.ambient = visual.vector(0.33, 0.33, 0.5)
    surfaceCanvas.up = visual.vector(0,0,1)
    surfaceCanvas.forward = visual.vector(0.33, -0.33, -0.15)
    surfaceCanvas.range = (tin.header.xMax - tin.header.xMin) * 0.55
    surfaceCanvas.caption = " *** COMMANDS ***\n\n 'r': run simulation \n 'p': pause "
    surfaceCanvas.caption += "\n 'u': move up (soil layer) \n 'd': move down (soil layer) "
    surfaceCanvas.caption += "\n 'l': load state \n 's': save state \n 'c': colorscale range"
    visual.label(canvas = surfaceCanvas, height = h, pos=visual.vector(cX, cY, Zlabel), text = "Surface")
    
    drawSurface(True)
    updateInterface()
    interface.bind('keydown', keyInput)


def drawColorScale():
    if C3DParameters.computeOnlySurface: 
        return
    
    step = (degreeMaximum - degreeMinimum) / nrColorLevels
    for i in range (nrColorLevels+1):
        degree = degreeMinimum + step * i
        c = getSEColor(degree, degreeMinimum, degreeMaximum)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = format(degree,".2f")


def updateColorScale():
    global degreeMinimum, degreeMaximum
    if C3DParameters.computeOnlySurface: 
        return
    
    rangeOk = False
    while not rangeOk:
        degreeMinimum = float(input("\nSet min. value (Degree of saturation):"))
        degreeMaximum = float(input("Set max. value (Degree of saturation):"))
        if degreeMaximum <= degreeMinimum:
            print("Wrong range!")
        else:
            rangeOk = True
            
    drawColorScale();
    drawSubSurface(False)


def updateLayer(s):
    global visualizedLayer
    if C3DParameters.computeOnlySurface: 
        return
    
    if s == 'd':
        if (visualizedLayer < C3DStructure.nrLayers-1):
            visualizedLayer += 1
    elif s == 'u':
        if (visualizedLayer > 1):
            visualizedLayer -= 1
             
    updateInterface()
    drawSubSurface(False)
          
                
def keyInput(evt):
    global isPause
    s = evt.key
    if s == 'r':
        isPause = False
    elif s == 'p':
        isPause = True
        print ("Pause...")
    elif s == 'd' or s == 'u':
        updateLayer(s)
    elif s == 's':
        isPause = True
        print ("Save State...")
        saveState()
    elif s == 'l':
        isPause = True
        print ("Load State...")
        if loadState(""):
            balance.initializeBalance()
            redraw()
    elif s == "c":
        if (isPause):
            updateColorScale()


def getNewTriangle(myColor, myCanvas, v):
    vert = []
    for i in range(3):
        vert.append(visual.vertex(pos = visual.vector(v[i,0], v[i,1], v[i,2] * tin.header.magnify)))
        vert[i].color = myColor
                     
    newTriangle = visual.triangle(canvas = myCanvas, vs=[vert[0],vert[1],vert[2]])
    return newTriangle


def drawSurface(isFirst):
    global surfaceCanves, surfaceTriangles
    maximum = 0.05
    for i in range(C3DStructure.nrTriangles):
        z = tin.C3DTIN[i].centroid[2]
        TINColor = getTINColor(z, tin.header)
        waterHeight = max(C3DCells[i].H - C3DCells[i].z, 0.0)
        waterColor = getSurfaceWaterColor(waterHeight, maximum)
        
        if waterHeight > maximum: 
            a = 1
        elif waterHeight < EPSILON_METER: 
            a = 0.0
        else: 
            a = max(0.2, waterHeight / maximum)
            
        c = a*waterColor + (1-a)*TINColor
        myColor = visual.vector(c[0], c[1], c[2])
        
        if (isFirst):
            newTriangle = getNewTriangle(myColor, surfaceCanvas, tin.C3DTIN[i].v)
            surfaceTriangles.append(newTriangle)
        else:
            surfaceTriangles[i].v0.color = myColor 
            surfaceTriangles[i].v1.color = myColor
            surfaceTriangles[i].v2.color = myColor 
 
 
def drawSubSurface(isFirst):
    global subSurfaceTriangles
    if C3DParameters.computeOnlySurface: 
        return
    
    depth = soil.depth[visualizedLayer] * 100
    layerLabel.text = "Degree of saturation " + format(depth,".1f")+"cm"
        
    for i in range(C3DStructure.nrTriangles):
        index = visualizedLayer * C3DStructure.nrTriangles + i
        c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
        myColor = visual.vector(c[0], c[1], c[2])
        
        if (isFirst):
            newTriangle = getNewTriangle(myColor, soilCanvas, tin.C3DTIN[i].v)
            subSurfaceTriangles.append(newTriangle)
        else:
            subSurfaceTriangles[i].v0.color = myColor
            subSurfaceTriangles[i].v1.color = myColor
            subSurfaceTriangles[i].v2.color = myColor    
    
    
def updateInterface():       
    timeLabel.text = "Time: " + str(int(balance.totalTime)) + " [s]"
    precLabel.text = "Prec: " + str(balance.currentPrec) + " [mm/hour]"
    storage = balance.currentStep.waterStorage
    flow = balance.currentStep.waterFlow
    timeStep = C3DParameters.currentDeltaT
    totalFlow = balance.allSimulation.waterFlow
    totalError = balance.allSimulation.MBE
    stepLabel.text = "Time step: " + str(timeStep) +" [s]"
    storageLabel.text = "Storage: " + format(storage,".2f") +" [m3]"
    flowLabel.text = "Flow: " + format(flow,".3f") + " [m3]"
    totalFlowLabel.text = "Total flow: " + format(totalFlow,".3f") + " [m3]"
    totalErrorLabel.text = "Total error: " + format(totalError,".3f") + " [m3]"
    
    
def redraw():
    updateInterface()
    drawSurface(False)
    if not C3DParameters.computeOnlySurface:
        drawSubSurface(False)
      
    