#PSP_GaussElimination.py
from __future__ import print_function, division
from PSP_readDataFile import readDataFile
from numpy import *

def GaussElimination(A,b):
    n=len(b)
    for m in range(0,n-1):
        for i in range(m+1,n):
            dummy=A[i,m]/A[m,m]
            A[i,m]=0.
            for j in range(m+1,n):
                A[i,j]-=(dummy*A[m,j])
            b[i]-=(dummy*b[m])
    # backward insertion
    x=zeros(n,float)
    x[n-1]=b[n-1]/A[n-1,n-1]
    for i in range(n-2,-1,-1):
        x[i]=b[i]
        for j in range(i+1,n):
            x[i]-=(A[i,j]*x[j])
        x[i]/=A[i,i]
    return (x)

def main():
    A, isFileOk = readDataFile("A.txt", 0, ',', False)
    if (isFileOk == False): print ("Incorrect format in row: ", A)
    x, isFileOk = readDataFile("b.txt", 0, ',', False)
    b = x[0,:]
    print("matrix A:\n",A)
    print("vector b:\n",b)
    
    x= GaussElimination(A,b)
    print("unknown vector x:\n",x) 

main()
