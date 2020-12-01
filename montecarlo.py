from math import *
import odlib
import numpy as np
import matplotlib.pyplot as plt

useLightspeedCorrection = True
radecAsDecimal = True

# Test Case:
# realElements = [1.056800391773215, .3442329516328222, 25.15526140011037,
#         236.2379969538250, 255.5046626217337, 139.5152750732308]
# t1 = 2458303.5
# t2 = 2458313.5
# t3 = 2458323.5
#
# realR1 = np.array([2.809885947589762E-01, -1.244695239934626E+00, 4.345462810378481E-01])
# realR2 = np.array([3.970631619077938E-01, -1.225073762366858E+00, 4.747424470178567E-01])
# realR3 = np.array([5.086027140699297E-01, -1.191426064461829E+00, 5.095072621707616E-01])
#
#
# ra1 = odlib.HMStoDeg(18, 41, 48.08)
# ra2 = odlib.HMStoDeg(18, 14, 30.53)
# ra3 = odlib.HMStoDeg(17, 54, 29.36)
# dec1 = odlib.DMStoDeg(36, 15, 14)
# dec2 = odlib.DMStoDeg(36, 9, 46)
# dec3 = odlib.DMStoDeg(34, 29, 32.2)
#
# R = np.array([[-2.068851717348834E-01, 9.132770535148356E-01, 3.958800809868246E-01],
#     [-3.688942164347761E-01, 8.690938354463554E-01, 3.767267647401627E-01],
#     [-5.204759869724843E-01, 8.004004162505617E-01, 3.469487530293208E-01]])
#
# eclipticObliquity = 23.4368815858

# Speed of light in AU / day
c = 173.145

realElements = [3.341858158593617E+00,6.312370547375413E-01,1.254126493887905E+01,
2.324514097688120E+02,5.116168365065160E+01,3.582977103336309E+02]

t1 = 2458667.6569132637
t2 = 2458671.6550316783
t3 = 2458677.731183194

realR1 = np.array([2.809885947589762E-01, -1.244695239934626E+00, 4.345462810378481E-01])
realR2 = np.array([3.970631619077938E-01, -1.225073762366858E+00, 4.747424470178567E-01])
realR3 = np.array([5.086027140699297E-01, -1.191426064461829E+00, 5.095072621707616E-01])


ra1 = 15.674761
dec1 = 9.704156
ra2 = 15.835688
dec2 = 11.152362
ra3 = 16.105651
dec3 = 12.936407

R = np.array([[-1.887642493130150E-01, 9.166674144548658E-01, 3.973362003879022E-01],
    [-2.547806549406440E-01, 9.031163132837493E-01, 3.914640974474277E-01],
    [-3.527852861550124E-01, 8.747684220168056E-01, 3.791725558768219E-01]])

eclipticObliquity = 23.4367541414

# Gets rho-hat from the RA and DEC
# Only omit last four args if the first two are RA in decimalized hours and
# Dec in decimalized degrees
def getRhoHat(rah, ram, ras=0, dech=0, decm=0, decs=0):
    # First two elements will be RA (hr) and DEC (deg) if radecAsDecimal
    if radecAsDecimal:
        raRad = odlib.degreesToRadians(rah*15)
        decRad = odlib.degreesToRadians(ram)
    else:
        raRad = odlib.degreesToRadians(odlib.HMStoDeg(rah, ram, ras))
        decRad = odlib.degreesToRadians(odlib.DMStoDeg(dech, decm, decs))
    rhox = cos(raRad) * cos(decRad)
    rhoy = sin(raRad) * cos(decRad)
    rhoz = sin(decRad)

    return np.array([rhox, rhoy, rhoz])

# Calculates the orbital elements of an asteroid given its position and
# velocity vectors
# Returns a, e, i, omega, w, and M
def orbitalElements(r, rDot):
    rMag = odlib.mag(r)
    vSqr = odlib.dot(rDot, rDot)
    a = 1 / (2 / rMag - vSqr)
    h = odlib.cross(r, rDot)
    hMag = odlib.mag(h)
    e = sqrt(1 - hMag**2 / a)
    i = acos(h[2] / hMag)
    cosOmega = -h[1] / (hMag * sin(i))
    sinOmega = h[0] / (hMag * sin(i))
    omega = odlib.sinCosToAngle(sinOmega, cosOmega)
    sinNu = (a * (1 - e**2) / hMag * (odlib.dot(r, rDot) / rMag)) / e
    cosNu = (a * (1 - e**2) / rMag - 1) / e
    nu = odlib.sinCosToAngle(sinNu, cosNu)
    cosU = (r[0] * cosOmega + r[1] * sinOmega) / rMag
    sinU = r[2] / (rMag * sin(i))
    u = odlib.sinCosToAngle(sinU, cosU)
    lowerOmega = u - nu
    E = acos((1 - rMag / a) / e)
    M = E - e * sin(E)

    return [a, e, i, omega, lowerOmega, M]

def od(t1, t2, t3, ra1, dec1, ra2, dec2, ra3, dec3, R):
    tau3 = odlib.timeToGaussian(t3, t2)
    tau1 = odlib.timeToGaussian(t1, t2)
    tau0 = tau3 - tau1

    realR2Dot = (realR3 - realR1) / tau0

    if radecAsDecimal:
        rhoHat1 = getRhoHat(ra1, dec1)
        rhoHat2 = getRhoHat(ra2, dec2)
        rhoHat3 = getRhoHat(ra3, dec3)
    else:
        rhoHat1 = getRhoHat(*ra1, *dec1)
        rhoHat2 = getRhoHat(*ra2, *dec2)
        rhoHat3 = getRhoHat(*ra3, *dec3)

    a = [tau3 / tau0, -1, -tau1 / tau0]

    rhoHats = np.array([rhoHat1, rhoHat2, rhoHat3])
    rho = rhoHats.copy()
    D0 = odlib.tri_product(rhoHat1, rhoHat2, rhoHat3)
    assert(D0 != 0)

    D1 = odlib.tri_product(rhoHat3, R[0], rhoHat2)
    D2 = odlib.tri_product(rhoHat3, R[1], rhoHat2)
    D3 = odlib.tri_product(rhoHat3, R[2], rhoHat2)
    rho[0] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[0] * D0)

    D1 = odlib.tri_product(rhoHat3, rhoHat1, R[0])
    D2 = odlib.tri_product(rhoHat3, rhoHat1, R[1])
    D3 = odlib.tri_product(rhoHat3, rhoHat1, R[2])
    rho[1] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[1] * D0)

    D1 = odlib.tri_product(rhoHat1, rhoHat2, R[0])
    D2 = odlib.tri_product(rhoHat1, rhoHat2, R[1])
    D3 = odlib.tri_product(rhoHat1, rhoHat2, R[2])
    rho[2] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[2] * D0)
    r = rho - R
    r2Dot = (r[2] - r[0]) / tau0

    tauAt1 = odlib.timeToGaussian(t1)
    tauAt2 = odlib.timeToGaussian(t2)
    tauAt3 = odlib.timeToGaussian(t3)

    def f(tau, tau2, r2, r2Dot):
        res = 1 - (tau - tau2)**2 / (2 * odlib.mag(r2)**3)
        res += odlib.dot(r2, r2Dot) / (2 * odlib.mag(r2)**5) * (tau - tau2)**3
        O4 = 3 * (odlib.dot(r2Dot, r2Dot) / odlib.mag(r2)**2 - 1 / odlib.mag(r2)**3)
        O4 += -15 * (odlib.dot(r2, r2Dot) / odlib.mag(r2)**2)**2 + 1 / odlib.mag(r2)**3
        O4 *= ((tau - tau2)**4 / (24 * odlib.mag(r2)**3))
        return res + O4
    def g(tau, tau2, r2, r2Dot):
        res = tau - tau2 - (tau - tau2)**3 / (6 * odlib.mag(r2)**3)
        res += (odlib.dot(r2, r2Dot) * (tau - tau2)**4) / (4 * odlib.mag(r2)**5)
        return res

    tiny = 1E-9
    maxIterations = 200
    iteration = 0
    # Minus 100 so the loop doesn't exit immediately
    r2Old = r[1] - 100
    r2DotOld = r2Dot
    while (np.any(abs(r[1] - r2Old) > tiny) or np.any(abs(r2Dot - r2DotOld) > tiny)) and iteration < maxIterations:
        r2Old = r[1].copy()
        r2DotOld = r2Dot.copy()

        if useLightspeedCorrection:
            # Lightspeed correction
            correctedT1 = t1 - odlib.mag(rho[0]) / c
            correctedT2 = t2 - odlib.mag(rho[1]) / c
            correctedT3 = t3 - odlib.mag(rho[2]) / c
            tauAt1 = odlib.timeToGaussian(correctedT1)
            tauAt2 = odlib.timeToGaussian(correctedT2)
            tauAt3 = odlib.timeToGaussian(correctedT3)

        rho = rhoHats.copy()
        D0 = odlib.tri_product(rhoHat1, rhoHat2, rhoHat3)
        assert(D0 != 0)

        D1 = odlib.tri_product(rhoHat3, R[0], rhoHat2)
        D2 = odlib.tri_product(rhoHat3, R[1], rhoHat2)
        D3 = odlib.tri_product(rhoHat3, R[2], rhoHat2)
        rho[0] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[0] * D0)

        D1 = odlib.tri_product(rhoHat3, rhoHat1, R[0])
        D2 = odlib.tri_product(rhoHat3, rhoHat1, R[1])
        D3 = odlib.tri_product(rhoHat3, rhoHat1, R[2])
        rho[1] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[1] * D0)

        D1 = odlib.tri_product(rhoHat1, rhoHat2, R[0])
        D2 = odlib.tri_product(rhoHat1, rhoHat2, R[1])
        D3 = odlib.tri_product(rhoHat1, rhoHat2, R[2])
        rho[2] *= (a[0] * D1 + a[1] * D2 + a[2] * D3) / (a[2] * D0)

        r = rho - R
        f1 = f(tauAt1, tauAt2, r[1], r2Dot)
        f3 = f(tauAt3, tauAt2, r[1], r2Dot)
        g1 = g(tauAt1, tauAt2, r[1], r2Dot)
        g3 = g(tauAt3, tauAt2, r[1], r2Dot)

        a[0] = g3 / (f1 * g3 - f3 * g1)
        a[2] = g1 / (f3 * g1 - f1 * g3)
        r[1] = a[0] * r[0] + a[2] * r[2]
        r2Dot = f3 / (g1 * f3 - g3 * f1) * r[0] + f1 / (g3 * f1 - g1 * f3) * r[2]
        iteration += 1

    assert(iteration < maxIterations)
    # print(r)
    r[1] = odlib.rotateVector(r[1], 0, odlib.degreesToRadians(eclipticObliquity), 0)
    r2Dot = odlib.rotateVector(r2Dot, 0, odlib.degreesToRadians(eclipticObliquity), 0)
    elements = orbitalElements(r[1], r2Dot)

    elements[2] = odlib.radiansToDegrees(elements[2])
    elements[3] = odlib.radiansToDegrees(elements[3])
    elements[4] = odlib.radiansToDegrees(elements[4])
    elements[5] = odlib.radiansToDegrees(elements[5])

    while elements[4] < 0:
        elements[4] += 360
    while elements[4] > 360:
        elements[4] -= 360

    return elements


ralist = declist = []
elementList = []
a_vals = []
e_vals = []
i_vals = []
omega_vals = []
w_vals = []
M_vals = []
for i in range(10000):
    randra1 = np.random.normal(ra1, 0.403/3600)
    randra2 = np.random.normal(ra2, 0.403/3600)
    randra3 = np.random.normal(ra3, 0.403/3600)
    randdec1 = np.random.normal(dec1, 0.441/3600)
    randdec2 = np.random.normal(dec2, 0.441/3600)
    randdec3 = np.random.normal(dec3, 0.441/3600)
    ralist.append([randra1, randra2, randra3])
    declist.append([randdec1, randdec2, randdec3])
    elementList.append(od(t1, t2, t3, randra1, randdec1, randra2, \
                            randdec2, randra3, randdec3, R))
    elements = elementList[i]
    a_vals.append(elements[0])
    e_vals.append(elements[1])
    i_vals.append(elements[2])
    omega_vals.append(elements[3])
    w_vals.append(elements[4])
    M_vals.append(elements[5])
    percentErrors = []
    for i in range(len(elements)):
        # I know this is over 80 characters
        # but I can't really make it shorter without shortening variable names
        percentErrors.append(100 * abs(elements[i] - realElements[i]) / realElements[i])

plt.hist(a_vals)
print(f"{odlib.mean(a_vals)} ± {odlib.stdev(a_vals)}")
plt.show()
plt.hist(e_vals)
print(f"{odlib.mean(e_vals)} ± {odlib.stdev(e_vals)}")
plt.show()
plt.hist(i_vals)
print(f"{odlib.mean(i_vals)} ± {odlib.stdev(i_vals)}")
plt.show()
plt.hist(omega_vals)
print(f"{odlib.mean(omega_vals)} ± {odlib.stdev(omega_vals)}")
plt.show()
plt.hist(w_vals)
print(f"{odlib.mean(w_vals)} ± {odlib.stdev(w_vals)}")
plt.show()
plt.hist(M_vals)
print(f"{odlib.mean(M_vals)} ± {odlib.stdev(M_vals)}")
plt.show()

    # print("Name\tExpected Value\t\tCalculated Value\tPercent Error")
    # print(f"a\t{realElements[0]}\t{elements[0]}\t{percentErrors[0]}%")
    # print(f"e\t{realElements[1]}\t{elements[1]}\t{percentErrors[1]}%")
    # print(f"i\t{realElements[2]}\t{elements[2]}\t{percentErrors[2]}%")
    # print(f"Omega\t{realElements[3]}\t{elements[3]}\t{percentErrors[3]}%")
    # print(f"w\t{realElements[4]}\t{elements[4]}\t{percentErrors[4]}%")
    # print(f"M\t{realElements[5]}\t{elements[5]}\t{percentErrors[5]}%")


# elements = od(t1, t2, t3, ra1, dec1, ra2, dec2, ra3, dec3, R)
# percentErrors = []
# for i in range(len(elements)):
#     # I know this is over 80 characters
#     # but I can't really make it shorter without shortening variable names
#     percentErrors.append(100 * abs(elements[i] - realElements[i]) / realElements[i])
#
# print("Name\tExpected Value\t\tCalculated Value\tPercent Error")
# print(f"a\t{realElements[0]}\t{elements[0]}\t{percentErrors[0]}%")
# print(f"e\t{realElements[1]}\t{elements[1]}\t{percentErrors[1]}%")
# print(f"i\t{realElements[2]}\t{elements[2]}\t{percentErrors[2]}%")
# print(f"Omega\t{realElements[3]}\t{elements[3]}\t{percentErrors[3]}%")
# print(f"w\t{realElements[4]}\t{elements[4]}\t{percentErrors[4]}%")
# print(f"M\t{realElements[5]}\t{elements[5]}\t{percentErrors[5]}%")
