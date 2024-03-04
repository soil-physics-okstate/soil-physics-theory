#PSP_triangulation.py
from PSP_gis import *
from PSP_triangulation import searchPosition

from re import split                                    #regular expression 
import sys
if sys.version_info < (3, 0):
    from Tkinter import *                           #2.7
    from tkFileDialog import askopenfilename            
else:
    from tkinter import *                           #3.x
    from tkinter.filedialog import askopenfilename
   
    
def openDTM():
    tkRoot = Tk()
    options = {}
    options['defaultextension'] = '.flt'
    options['filetypes'] = [('ESRI raster files', '.flt')]
    options['initialdir'] = 'data/'
    options['title'] = 'Open DTM'
    fileName = askopenfilename(**options)
    tkRoot.destroy()
    
    if (fileName != ""):
        header = loadEsriHeader(fileName[:-4] + ".hdr")
        dtm = loadEsriBinary(fileName, header)
        return(header, dtm)
    
    
def loadEsriHeader(fileName):
    header = Cheader()
    headerFile = open(fileName, "r")
    txtLine = headerFile.read()
    #separators: newline, one or more spaces, tab
    values = split('\n| +|\t', txtLine) 
    i = 0
    while (i < len(values)):
        tag = values[i].upper()
        if (tag == "NCOLS"):
            header.nrCols = int(values[i+1])
        elif (tag == "NROWS"):
            header.nrRows = int(values[i+1])
        elif (tag == "XLLCORNER"):
            header.xllCorner = float(values[i+1])
        elif (tag == "YLLCORNER"):
            header.yllCorner = float(values[i+1])
        elif (tag == "CELLSIZE"):
            header.cellSize = float(values[i+1]) 
        elif ((tag == "NODATA") or (tag == "NODATA_VALUE")):
            header.flag = float(values[i+1])
        i += 2
    return(header)
         
def loadEsriBinary(fileName, header):
    print ("load DTM data...")
    myFile = open(fileName, "rb")
    grid = np.fromfile(myFile, dtype=np.float32).reshape(
                header.nrCols, header.nrRows, order='F')
    computeMinMaxGrid(header, grid)
    return(grid) 


def isAdjacentIndex(t1, t2):
    shareVertices = 0
    for i in range(3):
        for j in range(3):
            if (t1[i] == t2[j]):
                shareVertices += 1
                if shareVertices == 2:
                    return True
    return False
    
def getNeighbours(triangleList): 
    nrTriangles = len(triangleList)  
    neighbourList = np.zeros((nrTriangles, 3), int)
    nrNeighbours = np.zeros(nrTriangles, int)
    for i in range(nrTriangles):
        index = 0
        neighbourList[i] = [NOLINK, NOLINK, NOLINK]
        j = 0
        while (j < nrTriangles) and (index < 3):
            if (nrNeighbours[j]<3):
                if isAdjacentIndex(triangleList[i], triangleList[j]):
                    if (j != i):
                        neighbourList[i, index] = j
                        nrNeighbours[j] += 1
                        index += 1
            j += 1
    return(neighbourList) 


def writeTIN(triangleList, pointList, header, dtm, 
             fnVertices, fnTriangles, fnNeighbours):
    print("Save TIN...")
    #save vertices
    myFile = open(fnVertices, "w")
    for i in range(len(pointList)):
        p = pointList[i]
        myFile.write(str(p[0]) +","+ str(p[1]) +","+ format(p[2],".1f") +"\n")
        
    #save triangle vertices
    triangleVertexList = np.zeros((len(triangleList),3), int)
    for i in range(len(triangleList)):
        for j in range(3):
            x = triangleList[i].v[j][0]
            pointIndex = searchPosition(x, pointList, 0, len(pointList)-1)
            while not(np.all(triangleList[i].v[j] == pointList[pointIndex])):
                pointIndex += 1
            triangleVertexList[i,j] = pointIndex
            
    myFile = open(fnTriangles, "w")
    for i in range(len(triangleVertexList)):
        t = triangleVertexList[i]
        myFile.write(str(t[0]) +","+ str(t[1]) +","+ str(t[2]) +"\n")
    
    neighbourList =  getNeighbours(triangleVertexList)
    myFile = open(fnNeighbours, "w")
    for i in range(len(neighbourList)):
        n = neighbourList[i]
        myFile.write(str(n[0]) +","+ str(n[1]) +","+ str(n[2]) +"\n")
            