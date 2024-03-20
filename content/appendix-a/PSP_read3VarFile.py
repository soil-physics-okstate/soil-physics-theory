#PSP_read3VarFile.py
from csv import reader

def read3VarFile(myFile, nRowHeader, myDelimiter, x1, x2, x3, printOnScreen):

	myReader = reader(open(myFile, "rt"), delimiter=myDelimiter)
	nRow = 0
	for row in myReader:
		nRow = nRow + 1
		# Header
		if (nRow <= nRowHeader):
			if (printOnScreen): print (row[0], row[1], row[2])
		# Values
		else:
			i = nRow - nRowHeader - 1
			x1.append(float(row[0]))
			x2.append(float(row[1]))
			x3.append(float(row[2]))
			if (printOnScreen):  print (x1[i], x2[i], x3[i])     