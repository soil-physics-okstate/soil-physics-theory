#PSP_CelsiusFahrenheit.py
from __future__ import print_function, division
import sys
if sys.version_info >= (3, 0):
    get_input = input
else:
    get_input = raw_input

def print_options():
    print ("Options:")
    print (" 'p' print options")
    print (" 'c' convert from Celsius")
    print (" 'f' convert from Fahrenheit")
    print (" 'q' quit the program")
 
def celsius_to_fahrenheit(c_temp):
    return (9.0 / 5.0) * c_temp + 32.0
 
def fahrenheit_to_celsius(f_temp):
    return (f_temp - 32.0) * 5.0 / 9.0
 
choice = "p"
while choice != "q":
    if choice == "c":
        temp = float(input("Celsius temperature: "))
        print ("Fahrenheit:", celsius_to_fahrenheit(temp), "\n")
    elif choice == "f":
        temp = float(input("Fahrenheit temperature: "))
        print ("Celsius:", fahrenheit_to_celsius(temp), "\n")
    print_options()
    choice = get_input()


       

