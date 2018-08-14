#hi from the outside
from numpy import *
from scipy import *
from pylab import *
import math
bloch = loadtxt("pbg3.txt", usecols=(0,), skiprows= 1, unpack =True)
freq = loadtxt("pbg3.txt", usecols=(1,), skiprows= 1, unpack =True)
symm = ['M',  r'$\Gamma$', 'X']
coord = [0, 10, 20, 34]
a =1
c = 3e8
const = a/c
ylabel('fa/c')
#ylim(0,1500)
#xlim(0,34)
axvline(x =10, color = 'black')
axvline(x =20, color = 'black')
xticks(coord, symm)
#title('Square Lattice PBD for EpsR=11.56, a=10e-6, r=0.18a')
plot(bloch, freq*a, 'ro')
show()

