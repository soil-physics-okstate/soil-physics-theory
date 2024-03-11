#PSP_triangle.py
from __future__ import division
from PSP_gis import *
from math import pi, acos

MINIMUMANGLE = pi/180.

class Crectangle:
    x0 = x1 = NODATA
    y0 = y1 = NODATA 
    
class Ccircle:
    def __init__(self, xc, yc, radiusSquared, isCorrect):
        if (isCorrect):
            self.x = xc
            self.y = yc 
            self.radiusSquared = radiusSquared
            self.radius = sqrt(radiusSquared)
        self.isCorrect =  isCorrect
       
class Cplane:
    # plane equation: ax + by + cz + d = 0
    a = b = c = d = NODATA
    
class Ctriangle:
    def __init__(self, v = np.zeros((3, 3), float)):
        self.v = v
        if (not np.all(v == 0.0)): 
            #centroid
            self.x = sum(v[:,0]) / 3.0
            self.y = sum(v[:,1]) / 3.0
            self.circle = getCircumCircle(v)
            self.isRefinedZ = False
            self.isRefinedAngle = False
    
def getRectangle(v):
    rect = Crectangle()
    x = v[:,0]
    y = v[:,1]
    rect.x0 = min(x[0], x[1], x[2])
    rect.y0 = min(y[0], y[1], y[2])
    rect.x1 = max(x[0], x[1], x[2])
    rect.y1 = max(y[0], y[1], y[2])
    return(rect)  

def getCircumCircle(v):
    if (getMinAngle(v) < MINIMUMANGLE):
        return Ccircle(NODATA, NODATA, NODATA, False)
    x = v[:,0]
    y = v[:,1]
    a = x[0]*x[0] + y[0]*y[0]
    b = x[1]*x[1] + y[1]*y[1]
    c = x[2]*x[2] + y[2]*y[2]
    d = 2*(x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]))
    xc = (a*(y[1] - y[2]) + b*(y[2] - y[0]) + c*(y[0] - y[1])) / d
    yc = (a*(x[2] - x[1]) + b*(x[0] - x[2]) + c*(x[1] - x[0])) / d
    dx = fabs(xc - x[0])
    dy = fabs(yc - y[0])
    radiusSquared = dx * dx + dy * dy
    return Ccircle(xc, yc, radiusSquared, True)
  
def getPlane(v):
    plane = Cplane()
    x = v[:,0]
    y = v[:,1]
    z = v[:,2]
    plane.a = (y[1]-y[0])*(z[2]-z[0]) - (y[2]-y[0])*(z[1]-z[0])
    plane.b = (x[2]-x[0])*(z[1]-z[0]) - (x[1]-x[0])*(z[2]-z[0])
    plane.c = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])
    plane.d = -(plane.a*x[0] + plane.b*y[0] + plane.c*z[0])
    return(plane)

def getZplane(plane, x, y):
    z = -(plane.a*x + plane.b*y + plane.d) / plane.c
    return(z)

# checks if a point is in circumcircle of the triangle.
def isInCircumCircle(point, triangle):
    dx = point[0] - triangle.circle.x
    dy = point[1] - triangle.circle.y
    if (((dx * dx) + (dy * dy)) <= triangle.circle.radiusSquared): return(True)
    else: return(False)

def sign(p, A, B):
    crossProductXY = (p[0]-A[0])*(B[1]-A[1]) - (p[1]-A[1])*(B[0]-A[0])
    if crossProductXY > 0.0: 
        return 1.0
    else:
        return -1.0
    
def isPointInside(point, v):
    s1 = sign(point, v[0], v[1])
    s2 = sign(point, v[1], v[2])
    s3 = sign(point, v[2], v[0])
    if (s1 == s2) and (s2 == s3): 
        return True
    else:
        return False
   
# Area = 1/2 |x1(y2-y3) - x2(y1-y3) + x3(y1 - y2)| 
def getArea2D(v):
    x = v[:,0]
    y = v[:,1]
    return 0.5 * fabs(x[0]*(y[1]-y[2]) - x[1]*(y[0]-y[2]) + x[2]*(y[0]-y[1]))

def magnitude(v):
    return(np.sqrt(v.dot(v)))

def getArea(v):
    return 0.5 * magnitude(np.cross(v[1] - v[0], v[2] - v[0]))

# Carnot equation
def getCosine(v0, v1, v2):
    a = distance2D(v0, v2)
    b = distance2D(v0, v1)
    c = distance2D(v2, v1)
    denom = 2.0 * b * c
    denom = max(denom, 0.00001)
    return (b*b + c*c - a*a) / denom

def getAngle(v0, v1, v2):
    cosine = getCosine(v0, v1, v2)
    cosine = min(cosine, 1.0)
    cosine = max(cosine, -1.0)
    return acos(cosine)  

def getPseudoAngle(v0, v1, v2):
    cosine = getCosine(v0, v1, v2)
    pseudoAngle = (cosine + 1.0)                #[ 0, 2]
    if sign(v0, v1, v2) > 0:  
        pseudoAngle *= -1                       #[-2, 2]
    return pseudoAngle

def getMinAngle(v):
    angle0 = getAngle(v[0],v[1],v[2])
    angle1 = getAngle(v[1],v[2],v[0])
    angle2 = getAngle(v[2],v[0],v[1])
    return min(angle0, angle1, angle2)      

#adjacent: 2 shared vertices with polygon v
def isAdjacent(u, v):
    nrSharedVertices = 0
    for i in range(3):
        for j in range(len(v)):
            if np.all(u[i] == v[j]):
                nrSharedVertices += 1
                if (nrSharedVertices == 2):
                    return True
    return False

def isCircumCircleLeft(triangle, point):
    xMax = triangle.circle.x + triangle.circle.radius
    currentX = point[0]
    return(currentX > xMax)

def hasVertexInDomain(triangle):
    for i in range(3):
        if (triangle.v[i][2] == NODATA): return(True)
    else: return(False)

def isInsideDTM(triangle, header, grid):
    zCircumCenter = getValueFromXY(header, grid, triangle.circle.x, triangle.circle.y)
    if (zCircumCenter != header.flag):
        return True
    else:
        return False
