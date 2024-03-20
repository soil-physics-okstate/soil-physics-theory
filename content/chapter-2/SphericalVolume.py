# PSP_basicProperties.py
from __future__ import print_function, division

def sphericalVolume(r):
        Volume = (4./3.)*3.14*(r**3)
        print ("\nThe volume of the sphere is = ", format(Volume,'.3f'))
        return

def cylinderVolume(r,h):
        cylinderVolume=3.14*(r**2)*h
        print ("\nThe volume of the cylinder is= ", format(cylinderVolume,'.3f'))

def main():
        r = float (input("input radius (m)= "))
        h=  float (input("input height (m)= "))
        sphericalVolume(r)
        cylinderVolume(r,h)
           
main()
