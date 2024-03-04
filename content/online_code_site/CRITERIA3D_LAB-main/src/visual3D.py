# visual3D.py
# ---------------------------------------------------------
# This module is part of the CRITERIA3D_LAB distribution
# https://github.com/ARPA-SIMC/CRITERIA3D_LAB
# ---------------------------------------------------------

import vpython as visual
from dataStructures import *
from color import *
from copy import copy
import waterBalance
import rectangularMesh
import soil
import crop
import importUtils

sliceRectangles = []
subSurfaceRectangles = []
visualizedLayer = 0
visualizedSlice = 0
nrColorLevels = 10
degreeMaximum = 1.0
degreeMinimum = 0.2
waterLevelMaximum = C3DParameters.pond
isWaterPotential = False
isRootVisualization = False
isPointVisualization = False
isPause = False
isComputeEquilibrium = False
colorScale = []


def initialize(totalWidth):
    global sliceCanvas, soilCanvas, interface
    global sliceLabel, layerLabel, timeLabel, precLabel, irrLabel, stepLabel, storageLabel
    global flowLabel, totalFlowLabel, totalErrorLabel
    global colorScale, visualizedSlice, visualizedLayer

    setAllColorScale()
    visualizedSlice = int(C3DStructure.nrRectanglesInYAxis / 2)

    # CANVAS DIMENSION
    interfaceWidth = int(totalWidth * 0.2)
    dx = int((totalWidth - interfaceWidth) / 2.0)
    dy = int(dx * 0.8)
    h = int(dy / 30)

    # CENTER
    cX = (rectangularMesh.header.xMin + rectangularMesh.header.xMax) * 0.5
    cY = (rectangularMesh.header.yMin + rectangularMesh.header.yMax) * 0.5
    cZ = rectangularMesh.header.zMax * rectangularMesh.header.magnify

    # RANGE
    rx = rectangularMesh.header.xMax - rectangularMesh.header.xMin
    ry = rectangularMesh.header.yMax - rectangularMesh.header.yMin
    rz = rectangularMesh.header.zMax - rectangularMesh.header.zMin + C3DStructure.gridDepth
    rangeX = rx * rectangularMesh.header.magnify
    rangeXY = max(rx, ry) * rectangularMesh.header.magnify
    rangeZ = rz * rectangularMesh.header.magnify

    # INTERFACE CANVAS
    interface = visual.canvas(width=interfaceWidth, height=dy, align="left")
    interface.background = visual.color.white
    interface.center = visual.vector(1, 0, 0)
    interface.range = 5

    timeLabel = visual.label(height=h, pos=visual.vector(0, 3, 0), text="", canvas=interface)
    stepLabel = visual.label(height=h, pos=visual.vector(0, 2, 0), text="", canvas=interface)
    precLabel = visual.label(height=h, pos=visual.vector(0, 1, 0), text="", canvas=interface)
    irrLabel = visual.label(height=h, pos=visual.vector(0, 0, 0), text="", canvas=interface)
    storageLabel = visual.label(height=h, pos=visual.vector(0, -1, 0), text="", canvas=interface)
    flowLabel = visual.label(height=h, pos=visual.vector(0, -2, 0), text="", canvas=interface)
    totalFlowLabel = visual.label(height=h, pos=visual.vector(0, -3, 0), text="", canvas=interface)
    totalErrorLabel = visual.label(height=h, pos=visual.vector(0, -4, 0), text="", canvas=interface)

    # COLOR LEGEND
    colorScale.clear()
    stepY = 10 / nrColorLevels
    for i in range(nrColorLevels + 1):
        label = visual.label(canvas=interface, pos=visual.vector(4, -5 + (i * stepY), 0), height=h,
                             background=visual.vec(0, 0, 0))
        colorScale.append(label)

    drawColorScale()

    # SURFACE CANVAS
    soilCanvas = visual.canvas(width=dx, height=dy, align="left")
    soilCanvas.background = visual.color.gray(1.0)
    soilCanvas.center = visual.vector(cX, cY, cZ)
    soilCanvas.ambient = visual.vector(0.2, 0.2, 0.2)
    soilCanvas.up = visual.vector(0, 0, 1)
    soilCanvas.forward = visual.vector(0, 1, -1)
    soilCanvas.range = rangeXY
    layerLabel = visual.label(color=visual.color.black, canvas=soilCanvas, height=h,
                              pos=visual.vector(cX, cY + rangeXY, cZ))

    drawSurface(True)

    # SLICE CANVAS
    sliceCanvas = visual.canvas(width=dx, height=dy, align="left")
    sliceCanvas.background = visual.color.gray(1.0)
    sliceCanvas.center = visual.vector(cX, cY, cZ - (rangeZ * 0.5))
    sliceCanvas.ambient = visual.vector(0.2, 0.2, 0.2)
    sliceCanvas.up = visual.vector(0, 0, 1)
    sliceCanvas.forward = visual.vector(0, 1, 0)
    sliceCanvas.range = max(rangeX, rangeZ)
    sliceLabel = visual.label(color=visual.color.black, canvas=sliceCanvas, height=h,
                              pos=visual.vector(cX, cY, cZ + (rangeZ * 0.2)))

    sliceCanvas.caption = " *** COMMANDS ***\n"
    sliceCanvas.caption += "\n '^': move up (soil layer) \n 'v': move down (soil layer) "
    sliceCanvas.caption += "\n '<': move left (soil slice) \n '>': move right (soil slice) "
    sliceCanvas.caption += "\n\n 'r': Run simulation \n 'p': Pause "
    sliceCanvas.caption += "\n 'e': compute hydrostatic Equilibrium"
    sliceCanvas.caption += "\n 's': Save state "
    sliceCanvas.caption += "\n 'l': Load state "
    sliceCanvas.caption += "\n "
    sliceCanvas.caption += "\n 'w': switch Water content - potential"
    sliceCanvas.caption += "\n 'd': view root Density"
    sliceCanvas.caption += "\n 'o': view Output points"
    sliceCanvas.caption += "\n 'c': set Colorscale range"

    drawSlice(True)
    updateInterface()
    interface.bind('keydown', keyInput)


def drawColorScale():
    if isWaterPotential:
        for i in range(nrColorLevels):
            colorScale[i].visible = True
        i = nrColorLevels
        c = getMatricPotentialColor(-0.01)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-0.1kPa"
        i -= 1
        c = getMatricPotentialColor(-0.1)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-1.0    "
        i -= 1
        c = getMatricPotentialColor(-0.5)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-5.0    "
        i -= 1
        c = getMatricPotentialColor(-1)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-10     "
        i -= 1
        c = getMatricPotentialColor(-2)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-20     "
        i -= 1
        c = getMatricPotentialColor(-3.3)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-33     "
        i -= 1
        c = getMatricPotentialColor(-5)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-50     "
        i -= 1
        c = getMatricPotentialColor(-10)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-100   "
        i -= 1
        c = getMatricPotentialColor(-30)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-300   "
        i -= 1
        c = getMatricPotentialColor(-100)
        colorScale[i].background = visual.vector(c[0], c[1], c[2])
        colorScale[i].text = "-1000 "
        for j in range(i):
            colorScale[j].visible = False
    elif isRootVisualization:
        for i in range(nrColorLevels):
            colorScale[i].visible = True
        step = 1 / nrColorLevels
        for i in range(nrColorLevels + 1):
            value = step * i
            c = getRootColor(value, 0, 1)
            colorScale[i].background = visual.vector(c[0], c[1], c[2])
            colorScale[i].text = format(value, ".2f")
    elif isPointVisualization:
        i = nrColorLevels
        colorScale[i].background = visual.color.white
        colorScale[i].text = "  Plant  "
        i -= 1
        colorScale[i].background = visual.color.cyan
        colorScale[i].text = "Dripper"
        i -= 1
        colorScale[i].background = visual.color.red
        colorScale[i].text = "Output "
        for j in range(i):
            colorScale[j].visible = False
    else:
        # water content
        for i in range(nrColorLevels):
            colorScale[i].visible = True
        step = (degreeMaximum - degreeMinimum) / nrColorLevels
        for i in range(nrColorLevels + 1):
            degree = degreeMinimum + step * i
            c = getSEColor(degree, degreeMinimum, degreeMaximum)
            colorScale[i].background = visual.vector(c[0], c[1], c[2])
            colorScale[i].text = format(degree, ".2f")


def updateColorScale():
    global degreeMinimum, degreeMaximum

    rangeOk = False
    while not rangeOk:
        degreeMinimum = float(input("\nSet min. value (Degree of saturation):"))
        degreeMaximum = float(input("Set max. value (Degree of saturation):"))
        if degreeMaximum <= degreeMinimum:
            print("Wrong range!")
        else:
            rangeOk = True

    drawColorScale()
    drawSurface(False)
    drawSlice(False)


def updateLayer(s):
    global visualizedLayer

    if s == 'down':
        if visualizedLayer < C3DStructure.nrLayers - 1:
            visualizedLayer += 1
    elif s == 'up':
        if visualizedLayer > 0:
            visualizedLayer -= 1

    updateInterface()
    drawSurface(False)


def updateSlice(s):
    global visualizedSlice

    if s == 'right':
        if visualizedSlice < C3DStructure.nrRectanglesInYAxis - 1:
            visualizedSlice += 1
    elif s == 'left':
        if visualizedSlice > 0:
            visualizedSlice -= 1

    updateInterface()
    drawSlice(False)


def keyInput(evt):
    global isPause, isComputeEquilibrium, isWaterPotential, isRootVisualization, isPointVisualization

    s = evt.key
    if s == 'r':
        isPause = False
    elif s == 'p':
        isPause = True
        print("Pause...")
    elif s == 'e':
        isComputeEquilibrium = True
    elif s == 'up' or s == 'down':
        updateLayer(s)
    elif s == 'left' or s == 'right':
        updateSlice(s)
    elif s == 's':
        isPause = True
        print("Save State...")
        fileName = importUtils.getStateFileName(True)
        if fileName != "":
            importUtils.saveCurrentModelState(fileName)
    elif s == 'l':
        isPause = True
        print("Load State...")
        fileName = importUtils.getStateFileName(False)
        if fileName != "":
            if importUtils.loadModelState(fileName):
                redraw()
            redraw()
    elif s == "c":
        if isPause:
            updateColorScale()
    elif s == "w":
        isWaterPotential = not isWaterPotential
        isRootVisualization = False
        isPointVisualization = False
        drawColorScale()
        redraw()
    elif s == "d":
        isRootVisualization = True
        isWaterPotential = False
        isPointVisualization = False
        drawColorScale()
        redraw()
    elif s == "o":
        isPointVisualization = True
        isWaterPotential = False
        isRootVisualization = False
        drawColorScale()
        redraw()


def getNewRectangle(myColor, myCanvas, v, computeNormal):
    vert = []
    for i in range(C3DStructure.nrVerticesPerRectangle):
        vert.append(visual.vertex(pos=visual.vector(v[i, 0], v[i, 1], v[i, 2] * rectangularMesh.header.magnify)))
        vert[i].color = myColor

    if computeNormal:
        normal = (vert[1].pos - vert[0].pos).cross(vert[2].pos - vert[1].pos).norm()
        for i in range(C3DStructure.nrVerticesPerRectangle):
            vert[i].normal = normal

    newRectangle = visual.quad(canvas=myCanvas, vs=[vert[0], vert[1], vert[2], vert[3]])
    return newRectangle


def drawSlice(isFirst):
    global sliceRectangles
    from crop import k_root, rootDensity, maxRootFactor
    from exportUtils import outputIndices

    firstIndex = visualizedSlice * C3DStructure.nrRectanglesInXAxis
    posY = C3DCells[firstIndex].y

    if isPointVisualization:
        var = "Output points"
    elif isRootVisualization:
        var = "Root density"
    elif isWaterPotential:
        var = "Water potential"
    else:
        var = "Degree of saturation"
    sliceLabel.text = var + " - slice at " + format(posY * 100, ".1f") + "cm"

    if isFirst:
        sliceRectangles.clear()

    for layer in range(C3DStructure.nrLayers):
        for x in range(C3DStructure.nrRectanglesInXAxis):
            index = firstIndex + x + (layer * C3DStructure.nrRectangles)
            i = layer * C3DStructure.nrRectanglesInXAxis + x
            if isRootVisualization:
                if index in outputIndices:
                    c = [1.0, 1.0, 1.0]
                else:
                    if layer == 0:
                        root = 0
                    else:
                        surfaceIndex = firstIndex + x
                        rootFactor = k_root[surfaceIndex] * rootDensity[surfaceIndex][layer] / soil.thickness[layer]
                        root = rootFactor / maxRootFactor
                    c = getRootColor(root, 0, 1)
            elif isPointVisualization:
                c = [0, 1, 0]
                if index in outputIndices:
                    c = [1, 0, 0]
            elif isWaterPotential:
                c = getMatricPotentialColor(C3DCells[index].H - C3DCells[index].z)
            else:
                c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)
            myColor = visual.vector(c[0], c[1], c[2])

            if isFirst:
                vertices = copy(rectangularMesh.C3DRM[firstIndex + x].v)
                for v in vertices[:2]:
                    v[2] = v[2] - soil.depth[layer] + (soil.thickness[layer] * 0.5)
                vertices[3] = vertices[0]
                vertices[2] = vertices[1]
                for v in vertices[2:]:
                    v[2] = v[2] - soil.thickness[layer]
                newRectangle = getNewRectangle(myColor, sliceCanvas, vertices, False)
                sliceRectangles.append(newRectangle)
            else:
                sliceRectangles[i].v0.color = myColor
                sliceRectangles[i].v1.color = myColor
                sliceRectangles[i].v2.color = myColor
                sliceRectangles[i].v3.color = myColor


def drawSurface(isFirst):
    global subSurfaceRectangles
    from exportUtils import outputSurfaceIndices, outputIndices

    maxWaterLevel = 0

    if isFirst:
        subSurfaceRectangles.clear()

    for i in range(C3DStructure.nrRectangles):
        index = visualizedLayer * C3DStructure.nrRectangles + i
        # color
        if visualizedLayer == 0:
            if isRootVisualization:
                c = getRootColor(crop.k_root[i], 0, 1.5)
            elif isPointVisualization:
                c = [0, 1, 0]
                if i in plantIndices:
                    c = [1, 1, 1]
                elif i in dripperIndices:
                    c = [0, 0.5, 1]
                elif i in outputSurfaceIndices:
                    c = [1, 0, 0]
            else:
                waterLevel = max(C3DCells[i].H - C3DCells[i].z, 0.0)
                maxWaterLevel = max(waterLevel, maxWaterLevel)
                c = getSurfaceWaterColor(waterLevel, waterLevelMaximum)
        else:
            if isWaterPotential:
                c = getMatricPotentialColor(C3DCells[index].H - C3DCells[index].z)
            elif isRootVisualization:
                if index in outputIndices:
                    c = [1.0, 1.0, 1.0]
                else:
                    rootFactor = crop.k_root[i] * crop.rootDensity[i][visualizedLayer] / soil.thickness[visualizedLayer]
                    root = rootFactor / crop.maxRootFactor
                    c = getRootColor(root, 0, 1)
            elif isPointVisualization:
                c = [0, 1, 0]
                if index in outputIndices:
                    c = [1, 0, 0]
            else:
                c = getSEColor(C3DCells[index].Se, degreeMinimum, degreeMaximum)

        myColor = visual.vector(c[0], c[1], c[2])

        if isFirst:
            newRectangle = getNewRectangle(myColor, soilCanvas, rectangularMesh.C3DRM[i].v, True)
            subSurfaceRectangles.append(newRectangle)
        else:
            subSurfaceRectangles[i].v0.color = myColor
            subSurfaceRectangles[i].v1.color = myColor
            subSurfaceRectangles[i].v2.color = myColor
            subSurfaceRectangles[i].v3.color = myColor

    # label
    if visualizedLayer == 0:
        if isPointVisualization:
            layerLabel.text = "Drippers and Output position"
        elif isRootVisualization:
            layerLabel.text = "Root factor"
        else:
            layerLabel.text = "Surface water level - max:" + format(maxWaterLevel * 1000, ".1f") + "mm"
    else:
        depth = soil.depth[visualizedLayer] * 100
        if isPointVisualization:
            var = "Output points"
        elif isRootVisualization:
            var = "Root density"
        elif isWaterPotential:
            var = "Water potential"
        else:
            var = "Degree of saturation"
        layerLabel.text = var + " - layer at " + format(depth, ".1f") + "cm"


def updateInterface():
    timeLabel.text = "Time: " + format(waterBalance.totalTime / 3600.0, ".3f") + " [h]"
    precLabel.text = "Rainfall: " + format(waterBalance.currentPrec, ".1f") + " [mm/hour]"
    irrLabel.text = "Irrigation: " + format(waterBalance.currentIrr, ".3f") + " [l/hour]"
    storage = waterBalance.currentStep.waterStorage
    flow = waterBalance.currentStep.waterFlow
    timeStep = C3DParameters.currentDeltaT
    totalFlow = waterBalance.allSimulation.waterFlow
    totalError = waterBalance.allSimulation.MBE
    stepLabel.text = "Time step: " + str(timeStep) + " [s]"
    storageLabel.text = "Storage: " + format(storage, ".5f") + " [m3]"
    flowLabel.text = "Current flow: " + format(flow * 1000, ".4f") + " [l]"
    totalFlowLabel.text = "Total flow: " + format(totalFlow * 1000, ".4f") + " [l]"
    totalErrorLabel.text = "Total error: " + format(totalError * 1000, ".4f") + " [l]"


def redraw():
    updateInterface()
    drawSlice(False)
    drawSurface(False)
