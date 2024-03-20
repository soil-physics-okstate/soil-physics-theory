from __future__ import print_function
from math import pi
from PSP_triangulation import triangulate
from PSP_utility import openDTM, writeTIN

import sys
if sys.version_info < (3, 0):
    import PSP_visual3D_27 as visual3D          #2.7
else:
    import PSP_visual3D as visual3D             #3.x
    

def main():
    #open DTM
    header, dtm = openDTM()
    print ("number of DTM points:", header.nrPoints)
    
    #show DTM
    visual3D.initialize(header)
    #visual3D.drawDTM(dtm, header)

    #set parameters
    compressionRate = 5
    angleThreshold = pi/3.0              # [rad] minimum angle for refinement
    initialPartition = False              # initial regular subdivision  
    
    #triangulation  
    pointList, triangleList = triangulate(header, dtm, compressionRate, angleThreshold, initialPartition)
    print ("number of triangles:", len(triangleList))
    visual3D.drawAllTriangles(triangleList, header)
    
    #show surface
    visual3D.drawAllSurfaces(triangleList, header)
    
    #write TIN
    writeTIN(triangleList, pointList, header, dtm, 
             "data/vertices.csv", "data/triangles.csv", "data/neighbours.csv")
    print("End.")      
main()
