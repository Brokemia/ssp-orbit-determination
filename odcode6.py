import odlib
import numpy as np
from math import sin, cos, sqrt, pi, asin

a = 3.333373705525415
e = 0.6300217569592357
i = odlib.degreesToRadians(12.607673026103942)
omega = odlib.degreesToRadians(232.6027678824531)
w = odlib.degreesToRadians(50.817782977450406)
M2 = odlib.degreesToRadians(1.6803544595919)
t2 = 2458671.6550316783

# Change these for different dates
t = 2458671.6550316783
eclipticObliquity = 23.4367540863
R = np.array([-2.547806549329457E-01, 9.031163132855038E-01, 3.914640974482125E-01])

# Finds the value of E in the equation M = E + e * sin(E)
# Uses Newton's method
def newtonFindE(e, realM):
    # An initial guess
    E = realM

    # 100 was chosen arbitrarily, but it works
    for i in range(100):
        M = E - e * sin(E)
        derivM = 1 - e * cos(E)
        E = E - (M - realM) / derivM
    return E

n = sqrt(1 / a**3)
T = t2 - M2 / (n * odlib.k)
M = M2 + n * odlib.k * (t - t2)
foundE = newtonFindE(e, M)

r = [a * cos(foundE) - a * e, a * sqrt(1 - e**2) * sin(foundE), 0]
r = odlib.rotateVectorZ(r, -w)
r = odlib.rotateVectorX(r, -i)
r = odlib.rotateVectorZ(r, -omega)
r = np.array(odlib.rotateVectorX(r, odlib.degreesToRadians(-eclipticObliquity)))
print(r)

rho = R + r
rhoHat = rho / odlib.mag(rho)

dec = asin(rhoHat[2])
cosRA = rhoHat[0] / cos(dec)
sinRA = rhoHat[1] / cos(dec)

RA = odlib.sinCosToAngle(sinRA, cosRA)

print("Predicted RA:\t[17, 42, 21.12]\nPredicted DEC:\t[31, 52, 26.7]")
print(f"Calculated RA: \t{odlib.RAdecimalToHMS(RA*180/pi)}\nCalculated DEC:\t{odlib.DECdecimalToDMS(dec*180/pi)}")
