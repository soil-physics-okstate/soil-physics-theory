#from __future__ import division
import numpy as np

NOLINK = -1

def linear(n, depth):
    z = np.zeros(n+2, float)
    dz = depth / n
    z[0] = 0
    z[1] = 0
    for i in range(1, n+1):
        z[i + 1] = z[i] + dz
    return z

def geometric(n, depth):
    z = np.zeros(n+2, float)  
    mySum = 0.0
    for i in range(1, n+1):
        mySum = mySum + i * i
    dz = depth / mySum
    z[0] = 0.0
    z[1] = 0.0
    for i in range(1, n+1):
        z[i + 1] = z[i] + dz * i * i
    return z


def getLayerIndex(z, depth):
    for i in range(1, len(z)-1):
        if (depth >= z[i]) and (depth <= z[i+1]):
            return(i)
    return(NOLINK)
    