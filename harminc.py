#Lbraries *******************************************************************
import numpy as np
from scipy import *
from pylab import *
from cmath import *
from numpy import ctypeslib
from ctypes import *
import harminv
for z in range (0,68):
    E1 = "EmptyLat%d" %z
    E = E1 + '.txt'
    print E
    i = 1
    P = np.zeros((5,20000), dtype=np.complex64)
    fmin = 0;
    fmax = 0.015
    dt = 0.01
    print i
    while (i < 11):
        print i
        E0R = loadtxt(E, usecols=(i,), skiprows= 1, unpack =True)
        i = i+1
        #print i
        E0I = loadtxt(E, usecols=(i,), skiprows= 1, unpack =True)
        l = i/2 - 1
        E0  = np.zeros(len(E0R), dtype=np.complex64)
        for m in range (0, len(E0R)):
            E0[m] = E0R[m] + 1j*E0I[m]
        P[l] = E0
        i = i+1
        #print i
    M = len(P[0])
    print E
    print M
    Mnew = M - 2000
    PSD = np.zeros((5, Mnew), dtype=np.complex64)
    for i in range (0,5):
        E0chop = P[i][2000: M]
        Mnew = len(E0chop)
        E0_inv = harminv.invert(E0chop, fmin = 0, fmax = 0.015, dt = 0.01)
        print E0_inv.freq
    #    PSD[i] = abs(E0_inv.freq)**2
    #    print M

    sumPSD  = np.zeros(Mnew)
    for j in range (0,Mnew):
        for i in range (0, 5):
            sumPSD[j] = sumPSD[j] + PSD[i][j]
            #plot(sumPSD)
            #show()
    dt = 1.66667e-11
    fs = 1/dt
    f = fs*arange(0,Mnew)/Mnew
    #plot(f, np.log(sumPSD))
    sender = np.log(sumPSD)
    #show()

    ff = open(E1 + "fft.txt", 'w')

    for k in range (0, Mnew):
            thang1 = str(f[k])
            thang2 = str(sender[k])
            ff.write(thang1)
            ff.write("\t")
            ff.write(thang2)
            ff.write("\n")
    ff.close()




