# Loads and calls a c library function
# Run using: python arrayuse.py

import numpy as np
from numpy import ctypeslib
from ctypes import *
import math
# Loads the c library
arraytest = cdll.LoadLibrary('./doWork.so.1.0')
# Provides a python function prototype for the c library function
arraytest.runsim.argtypes = [c_double, c_double, c_double, c_double, c_double]
# This is a wrapper function for the c library function
def runsim(k):
    #Call the library function
    arraytest.runsim(k[0], k[1], k[2], k[3], k[4])

# Create a kvector - actually nonsense
k = np.array([0.0, 0.0, 0.0, 50000, 0.0])
#i = np.array([0])
#Create an array of complex numbers - the dtype=np.complex64 is important


# Test the function
a =(0.5)
print a
GammaX = 2*math.pi/(a*math.sqrt(3))
#KxMax = 2*math.pi/(3*a)

print "Zone9"
print "X to Gamma"
k[0] = -GammaX/2
k[1] = -GammaX*math.sqrt(3)/2

step = GammaX/20
for i in range(0,20):
    runsim(k)
    k[0] = k[0] + step/2
    k[1] = k[1] + step*math.sqrt(3)/2
    k[4] +=1.0
    print k
    print "Zone9"




print "Gamma to J"
k[0] = 0
k[1] = 0
GammaJ = 2*GammaX/math.sqrt(3)
step = GammaJ/20
for i in range(0,20):
    runsim(k)
 #   k[0] = k[0] - step*math.sqrt(3)/2
    k[1] = k[1] - step
    k[4] +=1.0
    print k
    print "Zone9"

print "X to J"
XJ = GammaX/math.sqrt(3)
step = XJ/20

for i in range(0, 20):
    runsim(k)
    k[0] = k[0] - step*math.sqrt(3)/2
    k[1] = k[1] + step/2
    #k[0] = 2.0*math.pi/(math.sqrt(3.0)*a)*(i - 1)/20.0
    #k[1] = 2.0*math.pi/(3.0*a)*(i - 1)/20.0
    print k
    print "Zone9"
    k[4] +=1.0
