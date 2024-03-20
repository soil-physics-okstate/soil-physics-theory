#PSP_ThomasAlgorithm
from __future__ import print_function, division
from PSP_readDataFile import readDataFile
from PSP_ThomasAlgorithm import Thomas
import numpy as np

def main():
    myOutput, isFileOk = readDataFile("A.txt", 0, ',', False)
    if (isFileOk == False): 
        print ("Incorrect matrix format in row: ", myOutput)
        return (False)
    A = myOutput
    print("matrix A\n",A)
    
    myOutput, isFileOk = readDataFile("b.txt", 0, ',', False)
    d = myOutput[:,0]
    print("vector d\n",d)
    
    # initialize vectors
    n = len(d)   
    a = np.zeros(n, float)
    b = np.zeros(n, float)
    c = np.zeros(n, float)
    for i in range(n):
        if (i > 0): a[i] = A[i,i-1]
        b[i] = A[i,i]
        if (i < (n-1)): c[i] = A[i,i+1]
        
    x = Thomas(a, b, c, d)
    print("unknown vector x\n",x) 

main()
