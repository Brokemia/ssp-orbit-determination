from math import *
import copy as cp
import odlib
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as pat

# Calculates the determinant of a 2x2 matrix
def det2x2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]

# Calculates the determinant of a 3x3 matrix
def det3x3(m):
    d1 = m[0][0] * det2x2([[m[1][1], m[1][2]], [m[2][1], m[2][2]]])
    d2 = m[0][1] * det2x2([[m[1][0], m[1][2]], [m[2][0], m[2][2]]])
    d3 = m[0][2] * det2x2([[m[1][0], m[1][1]], [m[2][0], m[2][1]]])

    return d1 - d2 + d3

# Replaces column j in a 3x3 matrix with a 3x1 matrix
def spliceColumn(bigOlMatrix, smallOlMatrix, j):
    newMatrix = cp.deepcopy(bigOlMatrix)

    newMatrix[0][j] = smallOlMatrix[0][0]
    newMatrix[1][j] = smallOlMatrix[1][0]
    newMatrix[2][j] = smallOlMatrix[2][0]

    return newMatrix

# Uses Cramer's rule to solve a system of equations with 3 unknowns
def cramer(coefs, consts):
    a = det3x3(coefs)
    b1 = det3x3(spliceColumn(coefs, consts, 0))
    b2 = det3x3(spliceColumn(coefs, consts, 1))
    b3 = det3x3(spliceColumn(coefs, consts, 2))

    return [b1/a, b2/a, b3/a]

# Calculates the RA and DEC of an asteroid given the pixel location of the asteroid
# and the details of 6 stars (as Star objects)
# Also returns the uncertainty of the RA and DEC
def LSPR(x, y, stars):
    sumX = 0
    sumY = 0
    sumXY = 0
    sumXSqr = 0
    sumYSqr = 0
    sumRA = 0
    sumDEC = 0
    sumRAX = 0
    sumRAY = 0
    sumDECX = 0
    sumDECY = 0

    for star in stars:
        sumX += star.x
        sumY += star.y
        sumXY += star.x * star.y
        sumXSqr += star.x ** 2
        sumYSqr += star.y ** 2
        sumRA += star.RA
        sumDEC += star.DEC
        sumRAX += star.RA * star.x
        sumRAY += star.RA * star.y
        sumDECX += star.DEC * star.x
        sumDECY += star.DEC * star.y

    coefficients1 = [[len(stars),    sumX,    sumY],
                    [sumX, sumXSqr, sumXY],
                    [sumY, sumXY,   sumYSqr]]
    coefficients2 = [[len(stars),    sumX,    sumY],
                    [sumX, sumXSqr, sumXY],
                    [sumY, sumXY,   sumYSqr]]
    consts1 = [[sumRA], [sumRAX], [sumRAY]]
    consts2 = [[sumDEC], [sumDECX], [sumDECY]]

    # the best-fit plate coefficients
    resRA = cramer(coefficients1, consts1)
    resDEC = cramer(coefficients2, consts2)

    astRA = resRA[0] + resRA[1] * x + resRA[2] * y
    astDEC = resDEC[0] + resDEC[1] * x + resDEC[2] * y

    squaredRAResiduals = 0
    squaredDECResiduals = 0
    for star in stars:
        calcRA = resRA[0] + resRA[1] * star.x + resRA[2] * star.y
        calcDEC = resDEC[0] + resDEC[1] * star.x + resDEC[2] * star.y

        squaredRAResiduals += (star.RA - calcRA) ** 2
        squaredDECResiduals += (star.DEC - calcDEC) ** 2

    stdevRA = sqrt(squaredRAResiduals / (len(stars) - 3))
    stdevDEC = sqrt(squaredDECResiduals / (len(stars) - 3))

    return [astRA, astDEC, stdevRA, stdevDEC]


# A class used to hold details about stars on an image
class Star:
    def __init__(self, x, y, rah, ram, ras, decd, decam, decas, m):
        self.x = x
        self.y = y
        self.RAhr = rah
        self.RAmin = ram
        self.RAsec = ras
        self.RA = odlib.HMStoDeg(rah, ram, ras)
        self.DECdeg = decd
        self.DECamin = decam
        self.DECasec = decas
        self.DEC = odlib.DMStoDeg(decd, decam, decas)
        self.mag = m

# Reads the data necessary for a plate solve from a provided CSV file
def readPlateSolveDataFromCSV(path):
    f = open(path, 'r')
    csv_reader = csv.reader(f, delimiter=',')

    dataMatrix = []
    for row in csv_reader:
        dataMatrix.append([float(i) for i in row if i != ''])

    for i in range(len(dataMatrix)-1):
        dataMatrix[i] = Star(dataMatrix[i][0], dataMatrix[i][1], *odlib.RAdecimalToHMS(dataMatrix[i][2]), *odlib.DECdecimalToDMS(dataMatrix[i][3]), dataMatrix[i][4])

    return dataMatrix

csvdata = readPlateSolveDataFromCSV("platesolvedata.csv")
res = LSPR(csvdata[len(csvdata)-1][0], csvdata[len(csvdata)-1][1], csvdata[:len(csvdata)-1])
ra = odlib.RAdecimalToHMS(res[0])
print(f"Asteroid RA: {ra[0]:02d}:{ra[1]}:{ra[2]}")
dec = odlib.DECdecimalToDMS(res[1])
print(f"Asteroid Dec: {dec[0]}:{dec[1]}:{dec[2]}")
print("RA Uncertainty:", res[2] * 3600, "\nDec Uncertainty:", res[3] * 3600)

# Display graph of stars and asteroid
# dimmestStar = csvdata[0]
# for star in csvdata[:6]:
#     if star.mag < dimmestStar.mag:
#         dimmestStar = star
# for star in csvdata[:6]:
#     plt.plot(star.x, star.y, 'ro', markersize=(star.mag - dimmestStar.mag + 2)*4)
# ax = plt.gca()
# ax.add_patch(pat.Ellipse((csvdata[6][0], csvdata[6][1]), res[2]*200, res[3]*200, angle=0))
# plt.show()
