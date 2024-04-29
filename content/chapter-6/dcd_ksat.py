#PSP_Poisseulle.py
from visual import *
import numpy as np

scene.background=(1,1,1)
scene.title = ""

l = 5.     
R = 1.     
D_P = 40.   
eta = 1.   

angles=arange(1.1*pi,2.1*pi,pi/20.)
n = 500
for i in range(n):
 spring=curve(color=(0,1,1), radius=0.06)
 for phi in angles:
     spring.append(pos=(l*(float(i)/float(n)-0.55), R*cos(phi), R*sin(phi)))

for i in range(11):
  x = 0
  y = (float(i)/11.-0.5)*2.*R
  r = np.sqrt(x*x+y*y)
  if (r < R):
   arrow(pos=(0,y,x),axis = (D_P * (R*R-r*r)/(eta*l*4.),0,0),
                       shaftwidth = 0.03, color=(0,0,0))