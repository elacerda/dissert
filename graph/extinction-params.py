import sys
import os
import numpy as np 

# Cardelli et al. 1989, UV extinction, eqs. 4a and 4b.

def Fa(x):
    return -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3

def Fb(x):
    return 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3

# Eq. 4a
def a(x):
    return 1.752 - (0.316*x) - 0.104 / ((x-4.67)**2 + 0.341) + (Fa(x) if (x > 5.9) else 0.0) 

def b(x):
    return -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263) + (Fb(x) if (x > 5.9) else 0.0)


RV = 3.1

# Galex filters
FUV = 1528.0 # \AA
NUV = 2271.0 # \AA

xf = 10000.0 / FUV
xn = 10000.0 / NUV

print '\lambda_{FUV} = %f, \lambda_{NUV} = %f' % (FUV, NUV)
print 'x_{FUV} = %f, x_{NUF} = %f' % (xf, xn)

print 'A_{FUV} = %f E(B-V)' % (RV*a(xf) + b(xf))
print 'A_{NUV} = %f E(B-V)' % (RV*a(xn) + b(xn))

