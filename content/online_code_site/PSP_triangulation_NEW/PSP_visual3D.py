#PSP_visual3D.py
import vpython as visual
#from vpython import *
import numpy as np
from math import sqrt
from PSP_gis import getPointFromColRow
from PSP_color import setColorScaleDtm, getDTMColor

pointListDTM = []
triangleList3D = []

WINDOWS_SIZE = 600

def initialize(header):
    cX = header.xllCorner + header.nrCols * header.cellSize * 0.5
    cY = header.yllCorner + header.nrRows * header.cellSize * 0.5
    cZ = header.zMin + header.dz * 0.5
    
    visual.scene.width = WINDOWS_SIZE
    visual.scene.height = WINDOWS_SIZE
    visual.scene.background = visual.color.white
    visual.scene.ambient = visual.vector(0.33, 0.33, 0.5)
    visual.scene.center = visual.vector(cX, cY, cZ*header.magnify)
    visual.scene.up = visual.vector(0, 0, 1)
    visual.scene.forward = visual.vector(0, 1, -1)
    visual.scene.range = sqrt(header.nrPoints) * header.cellSize

    setColorScaleDtm()
    

def drawDTM(dtm, header):
    for col in range(header.nrCols):
        for row in range(header.nrRows):
            if dtm[col, row]!= header.flag:
                p = getPointFromColRow(header, dtm, col, row)
                c = getDTMColor(p[2], header)
                myColor = visual.vector(c[0],c[1],c[2])
                myPos = visual.vector(p[0], p[1], p[2] * header.magnify)
                myPoint = visual.points(pos = myPos, radius = 1, color=myColor) 
                pointListDTM.append(myPoint)
        
    print ("press any key to continue")
    visual.scene.waitfor('keydown')


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


def cleanAll():
    cleanAllPointsDTM()
    cleanAllTriangles()


def addTriangle(index, triangle, header):
    v = []
    myRadius = visual.scene.range / (WINDOWS_SIZE * 1.5)
    for i in range(3):
        v.append(visual.vector(triangle.v[i,0], triangle.v[i,1], triangle.v[i,2] * header.magnify))
    myTriangle = visual.curve(pos=[v[0], v[1] ,v[2], v[0]], color = visual.color.black, radius = myRadius)
    triangleList3D.insert(index, myTriangle)
    
         
def drawAllTriangles(triangleList, header):
    cleanAll()
    for i in range(len(triangleList)):
        addTriangle(i, triangleList[i], header)
        
    print ("press any key to continue")
    visual.scene.waitfor('keydown')


def drawSurface(v, header):
    z = sum(v[:,2]) / 3.0
    c = getDTMColor(z, header)
    myColor = visual.vector(c[0], c[1], c[2])
    vert = []
    for i in range(3):
        vert.append(visual.vertex(pos = visual.vector(v[i,0], v[i,1], v[i,2] * header.magnify), color = myColor))
        
    visual.triangle(vs=[vert[0], vert[1], vert[2]])
        
    
def drawAllSurfaces(triangleList, header):
    cleanAll()
    for i in range(len(triangleList)):
        drawSurface(triangleList[i].v, header)


