from math import *
import odlib
import numpy as np

# Gets rho-hat from the RA and DEC
def getRhoHat(rah, ram, ras, dech, decm, decs):
    raRad = odlib.degreesToRadians(odlib.HMStoDeg(rah, ram, ras))
    decRad = odlib.degreesToRadians(odlib.DMStoDeg(dech, decm, decs))
    rhox = cos(raRad) * cos(decRad)
    rhoy = sin(raRad) * cos(decRad)
    rhoz = sin(decRad)

    return np.array([rhox, rhoy, rhoz])

# R Vector:
# 2458303.500000000 = A.D. 2018-Jul-04 00:00:00.0000 TDB
#  X =-2.068851717348834E-01 Y = 9.132770535148356E-01 Z = 3.958800809868246E-01
#  VX=-1.654924029723420E-02 VY=-2.941750564327286E-03 VZ=-1.364862618215202E-03
#  LT= 5.871738278468436E-03 RG= 1.016659967384694E+00 RR= 1.936116458703009E-04
# 2458313.500000000 = A.D. 2018-Jul-14 00:00:00.0000 TDB
#  X =-3.688942164347761E-01 Y = 8.690938354463554E-01 Z = 3.767267647401627E-01
#  VX=-1.577976012898477E-02 VY=-5.468484327181214E-03 VZ=-2.458566706825568E-03
#  LT= 5.870978706896408E-03 RG= 1.016528451643865E+00 RR= 1.399157581197858E-04
# 2458323.500000000 = A.D. 2018-Jul-24 00:00:00.0000 TDB
#  X =-5.204759869724843E-01 Y = 8.004004162505617E-01 Z = 3.469487530293208E-01
#  VX=-1.454449550843527E-02 VY=-7.831479375401285E-03 VZ=-3.481186523108828E-03
#  LT= 5.866941738269600E-03 RG= 1.015829472193859E+00 RR= 9.248400657528764E-05

# Results:
# 2458303.500000000 = A.D. 2018-Jul-04 00:00:00.0000 TDB
#  X = 2.809885947589762E-01 Y =-1.244695239934626E+00 Z = 4.345462810378481E-01
#  VX= 1.179693344069725E-02 VY= 1.227514996916573E-03 VZ= 4.285330123314827E-03
#  LT= 7.785285815422619E-03 RG= 1.347980452775323E+00 RR= 2.707083724359027E-03
# 2458313.500000000 = A.D. 2018-Jul-14 00:00:00.0000 TDB
#  X = 3.970631619077938E-01 Y =-1.225073762366858E+00 Z = 4.747424470178567E-01
#  VX= 1.139883453711423E-02 VY= 2.679831160753917E-03 VZ= 3.750852820951060E-03
#  LT= 7.927086204450063E-03 RG= 1.372532429046544E+00 RR= 2.203048484238906E-03
# 2458323.500000000 = A.D. 2018-Jul-24 00:00:00.0000 TDB
#  X = 5.086027140699297E-01 Y =-1.191426064461829E+00 Z = 5.095072621707616E-01
#  VX= 1.089183150789963E-02 VY= 4.033874652042018E-03 VZ= 3.199511099610070E-03
#  LT= 8.039742324266541E-03 RG= 1.392038231530674E+00 RR= 1.698032246530710E-03

# Speed of light in AU / day
c = 173.145

t1 = 2458303.5
t2 = 2458313.5
t3 = 2458323.5

realR1 = np.array([3.959230509325300E-01, -1.225341042642005E+00, 4.743670942970302E-01])
#realR1 = np.array([3.959230509325300E-01, -1.312920820178346E+00, -5.218937513537718E-02])
realR2 = np.array([3.970631619077938E-01, -1.225073762366858E+00, 4.747424470178567E-01])
#realR3 = np.array([3.982028174780757E-01, -1.312727481764843E+00, -5.128791215021865E-02])
realR3 = np.array([3.982028174780757E-01, -1.224805076733292E+00, 4.751172648112766E-01])


ra1 = [18, 41, 48.08]
ra2 = [18, 14, 30.53]
ra3 = [17, 54, 29.36]
dec1 = [36, 15, 14]
dec2 = [36, 9, 46]
dec3 = [34, 29, 32.2]

R = np.array([[-2.068851717348834E-01, 9.132770535148356E-01, 3.958800809868246E-01],
    [-3.688942164347761E-01, 8.690938354463554E-01, 3.767267647401627E-01],
    [-5.204759869724843E-01, 8.004004162505617E-01, 3.469487530293208E-01]])

tau3 = odlib.timeToGaussian(t3, t2)
tau1 = odlib.timeToGaussian(t1, t2)
tau0 = tau3 - tau1

realR2Dot = (realR3 - realR1) / (odlib.timeToGaussian(2458313.6, t2) - odlib.timeToGaussian(2458313.4, t2))

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

r[1] = odlib.rotateVector(r[1], 0, odlib.degreesToRadians(23.4368815858), 0)
r2Dot = odlib.rotateVector(r2Dot, 0, odlib.degreesToRadians(23.4368815858), 0)

percentErrorR2 = abs(r[1] - realR2) / realR2 * 100
percentErrorR2Dot = abs(r2Dot - realR2Dot) / realR2Dot * 100
print('iterations', iteration)
print(f"Expected r: {realR2}\nCalulated r: {r[1]}")
print(f"Percent difference for each component: {percentErrorR2}\n")
print(f"Expected rdot: {realR2Dot}\nCalulated rdot: {r2Dot}")
print(f"Percent difference for each component: {percentErrorR2Dot}")
