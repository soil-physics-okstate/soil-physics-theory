#PSP_read6VarFile.py
from csv import reader

def read6VarFile(myFile, nRowHeader, myDelimiter, x1, x2, x3,x4,x5,x6, printOnScreen):

        myReader = reader(open(myFile, "rt"), delimiter=myDelimiter)
        nRow = 0
        for row in myReader:
            nRow = nRow + 1
            # Header
            if (nRow <= nRowHeader):
                if (printOnScreen): print (row[0], row[1], row[2],row[3],row[4],row[5])
            # Values
            else:
                i = nRow - nRowHeader - 1
                x1.append(float(row[0]))
                x2.append(float(row[1]))
                x3.append(float(row[2]))
                x4.append(float(row[3]))
                x5.append(float(row[4]))
                x6.append(float(row[5]))
                if (printOnScreen):  print (x1[i], x2[i], x3[i],x4[i],x5[i],x6[i])        
