# PSP_basicProperties.py
from __future__ import print_function, division

def sphericalVolume(r):
        Volume = (4./3.)*3.14*(r**3)
        print ("\nThe volume of the sphere is = ", format(Volume,'.6f'))
        return

def cylinderVolume(r,h):
        cylinderVolume=3.14*(r**2)*h
        print ("\nThe volume 1of the cylinder is= ", format(cylinderVolume,'.6f'))
        return

def massWater(r,h):
        massWater=3.14*(r**2)*h * 1000.
        print ("\nThe mass is (kg)= ", format(massWater,'.3f'))
        return
        
def main():
        r = float (input("input radius (m)= "))
        h=  float (input("input height (m)= "))
        sphericalVolume(r)
        cylinderVolume(r,h)
        massWater(r,h)
main()
