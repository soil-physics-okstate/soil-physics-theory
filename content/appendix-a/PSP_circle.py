# PSP_circle.py
from __future__ import print_function, division
from math import sqrt
import sys
if sys.version_info >= (3,0):
    from tkinter import *                           #3.x
else:
    from Tkinter import *                           #2.7

myCenterX = 150
myCenterY = 150
myRadius = 50
myRoot = Tk()
myText = StringVar()

def drawCircle(myCanvas,x,y,r):
    myCanvas.create_oval(x-r, y-r, x+r, y+r, width=1, fill='')

def ButtonLeftHandler(myEvent):
    dX = myEvent.x - myCenterX
    dY = myEvent.y - myCenterY
    distance = sqrt(dX*dX + dY*dY)
    if (distance <= myRadius):
        myText.set("Inside!")
    else:
        myText.set("Outside!")

def main():
    myCanvas = Canvas(width=300, height=300, bg='white')
    myCanvas.pack(expand=YES, fill=BOTH)
    myLabel = Label(myRoot, textvariable = myText, font=("Helvetica", 18))
    myLabel.pack()
    myText.set("Try to click the mouse button in the canvas...")
    myCircle = drawCircle(myCanvas, myCenterX, myCenterY, myRadius)
    myRoot.bind_class ("Canvas", "<Button-1>", ButtonLeftHandler)
    myRoot.mainloop()
main()
