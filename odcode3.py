from math import *
import odlib

# Gets rho-hat from the RA and DEC
def getRhoHat(rah, ram, ras, dech, decm, decs):
    raRad = odlib.degreesToRadians(odlib.HMStoDeg(rah, ram, ras))
    decRad = odlib.degreesToRadians(odlib.DMStoDeg(dech, decm, decs))
    rhox = cos(raRad) * cos(decRad)
    rhoy = sin(raRad) * cos(decRad)
    rhoz = sin(decRad)

    return [rhox, rhoy, rhoz]

# Gets rho-hat from rho
def rhoHatFromRho(rho):
    rhoMag = odlib.mag(rho)
    return [rho[0] / rhoMag, rho[1] / rhoMag, rho[2] / rhoMag]

rah = 18
ram = 41
ras = 48.08
dech = 36
decm = 15
decs = 14

realRho = [7.410342519834037E-02, -4.015610850300696E-01, 2.994571559854200E-01]
rho1 = [getRhoHat(rah, ram, ras, dech, decm, decs), rhoHatFromRho(realRho)]

rah = 18
ram = 14
ras = 30.53
dech = 36
decm = 9
decs = 46

realRho = [2.816894792575548E-02, -4.437310659481813E-01, 3.249880756106725E-01]
rho2 = [getRhoHat(rah, ram, ras, dech, decm, decs), rhoHatFromRho(realRho)]

rah = 17
ram = 54
ras = 29.36
dech = 34
decm = 29
decs = 32.2

realRho = [-1.187327019443038E-02, -4.953819751748095E-01, 3.404904429787050E-01]
rho3 = [getRhoHat(rah, ram, ras, dech, decm, decs), rhoHatFromRho(realRho)]

print(f"Rho-hat 1:\n\tFrom RA/DEC: {rho1[0]}\n\tFrom Rho: {rho1[1]}")
print(f"\tPercent Error X: {abs((rho1[0][0] - rho1[1][0]) / rho1[1][0]) * 100}%")
print(f"\tPercent Error Y: {abs((rho1[0][1] - rho1[1][1]) / rho1[1][1]) * 100}%")
print(f"\tPercent Error Z: {abs((rho1[0][2] - rho1[1][2]) / rho1[1][2]) * 100}%\n")
print(f"Rho-hat 2:\n\tFrom RA/DEC: {rho2[0]}\n\tFrom Rho: {rho2[1]}\n")
print(f"\tPercent Error X: {abs((rho2[0][0] - rho2[1][0]) / rho2[1][0]) * 100}%")
print(f"\tPercent Error Y: {abs((rho2[0][1] - rho2[1][1]) / rho2[1][1]) * 100}%")
print(f"\tPercent Error Z: {abs((rho2[0][2] - rho2[1][2]) / rho2[1][2]) * 100}%\n")
print(f"Rho-hat 3:\n\tFrom RA/DEC: {rho3[0]}\n\tFrom Rho: {rho3[1]}\n")
print(f"\tPercent Error X: {abs((rho3[0][0] - rho3[1][0]) / rho3[1][0]) * 100}%")
print(f"\tPercent Error Y: {abs((rho3[0][1] - rho3[1][1]) / rho3[1][1]) * 100}%")
print(f"\tPercent Error Z: {abs((rho3[0][2] - rho3[1][2]) / rho3[1][2]) * 100}%\n")
