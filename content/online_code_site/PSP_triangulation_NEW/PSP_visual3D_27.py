#PSP_visual3D.py
from __future__ import print_function, division
import numpy as np
from math import sqrt
from copy import copy
from PSP_gis import getPointFromColRow
from PSP_color import setColorScaleDtm, getDTMColor
import visual

pointListDTM = []
triangleList3D = []

def initialize(header):
    setColorScaleDtm()
    
    cX = header.xllCorner + header.nrCols * header.cellSize * 0.5
    cY = header.yllCorner + header.nrRows * header.cellSize * 0.5
    cZ = header.zMin + header.dz * 0.5
    
    global TINScene
    TINScene = visual.display(x = 0, y = 0, width = 600, height = 600)
    TINScene.title = "TIN"
    TINScene.background = visual.color.white
    TINScene.ambient = 0.33
    TINScene.center = (cX, cY, cZ*header.magnify)
    TINScene.up = (0,0,1)
    TINScene.forward = (0, 1, -1)
    TINScene.range = sqrt(header.nrPoints) * header.cellSize


def drawDTM(dtm, header):
    for col in range(header.nrCols):
        for row in range(header.nrRows):
            if dtm[col, row]!= header.flag:
                p = getPointFromColRow(header, dtm, col, row)
                myColor = getDTMColor(p[2], header)
                myPos = ([p[0], p[1], p[2] * header.magnify])
                myPoint = visual.points(display = TINScene, pos = myPos, size = 2.0, color=myColor) 
                pointListDTM.append(myPoint)
    print ("press a key to continue")
    TINScene.kb.getkey()


def cleanAllPointsDTM():      
    while (len(pointListDTM) > 0):
        pointListDTM[0].visible = False
        myPoint = pointListDTM.pop(0)
        del myPoint
 

def delTriangle(i):
    triangleList3D[i].visible = False 
    myTriangle = triangleList3D.pop(i)
    del myTriangle

     
def cleanAllTriangles():
    lastIndex = len(triangleList3D) - 1
    for i in np.arange(lastIndex, -1, -1):
        delTriangle(i)
 
 
def drawTriangle(index, triangle, header):
    v = np.zeros((3, 3), float)
    for i in range(3):
        v[i] = triangle.v[i]
    v[:,2] *= header.magnify
    myTriangle = visual.curve(pos=[v[0], v[1], v[2], v[0]], color=visual.color.black)
    triangleList3D.insert(index, myTriangle)  


def drawAllTriangles(triangleList, header):
    cleanAllPointsDTM()
    cleanAllTriangles()
    for i in range(len(triangleList)):
        drawTriangle(i, triangleList[i], header)
    print ("press a key to continue")
    TINScene.kb.getkey()
    
    
def drawSurface(myColor, v):
    mySurface = visual.faces(display = TINScene, color=myColor, pos = v)
    mySurface.make_twosided()
    mySurface.make_normals()
    mySurface.smooth()


def drawAllSurfaces(triangleList, header):
    cleanAllPointsDTM()
    cleanAllTriangles()
    nrTriangles = len(triangleList)
    for i in range(nrTriangles):
        v = copy(triangleList[i].v)
        z = sum(v[:,2]) / 3.0
        myColor = getDTMColor(z, header)
        v[:,2] *= header.magnify
        drawSurface(myColor, v)

