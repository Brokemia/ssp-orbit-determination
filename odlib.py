# Various utility functions to help with orbit determination
# Ricky Weerts
from math import *

# Calculates the dot product of two vectors
def dot(vector1, vector2):
    res = 0

    # For each component in the vectors, multiply them together
    # Then add them to the total stored in res
    for i in range(0, len(vector1)):
        res += vector1[i] * vector2[i]
    return res


# Calculates the cross product of two vectors
def cross(vector1, vector2):
    res_vector = [0,0,0]

    # Calculate each component
    res_vector[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1]
    res_vector[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2]
    res_vector[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0]

    return res_vector


# Calculates the triple product of three vectors
def tri_product(a, b, c):
    return dot(a, cross(b, c))

# Calculates the mean of a list of values
def mean(x):
    total = 0
    for i in x:
        total += i

    return total/len(x)

# Calculates the standard deviation of a list of values
def stdev(x):
    avg = mean(x)
    squaredDiffs = 0
    for i in x:
        squaredDiffs += (i - avg)**2

    return sqrt(squaredDiffs / (len(x)-1))

# Converts hours, minutes, and seconds to degrees for RA
def HMStoDeg(h, m, s):
    time = h + m/60 + s/3600
    return time / 24 * 360

# Converts degrees, arcminutes, and arcseconds to degrees for declination
def DMStoDeg(d, m, s):
    return d + copysign(1, d)*m/60 + copysign(1, d)*s/3600

# Converts RA to hours, minutes, and seconds and then prints the result
def RAdecimalToHMS(ra):
    ra = ra /360 *24
    h = floor(ra)
    m = (ra - h) * 60
    s = (m - floor(m)) * 60
    return [h, floor(m), s]

# Converts declination to degrees, arcminutes, and arcseconds
def DECdecimalToDMS(dec):
    d = int(dec)
    m = abs(dec - d) * 60
    s = abs((m - floor(m)) * 60)
    return [d, floor(m), s]

# Returns the magnitude of a vector
def mag(v):
    return sqrt(v[0]**2+v[1]**2+v[2]**2)

# Converts a sine and cosine of an angle to the angle in the correct quadrant
def sinCosToAngle(s, c):
    angle = asin(s)

    if c < 0:
        angle = pi - angle

    if angle < 0:
        angle = 2*pi + angle

    return angle

# Take two sides and one angle of a spherical triangle and return
# the other side and angles
def sphericalTriangleFromTwoSidesAndAngle(a, c, B):
    cosb = cos(a) * cos(c) + sin(a) * sin(c) * cos(B)
    A = acos((cos(a) - cosb * cos(c)) / (sin(acos(cosb)) * sin(c)))
    sinb = (cos(a) - cosb * cos(c)) / (sin(c) * cos(A))
    b = sinCosToAngle(sinb, cosb)
    C = acos((cos(c) - cos(a) * cosb) / (sin(a) * sin(b)))
    return [b, A, C]

# Converts degrees to radians
def degreesToRadians(d):
    return d * (pi / 180)

# Converts radians to degrees
def radiansToDegrees(r):
    return r * (180 / pi)

# Performs a single rotation when given a rotation matrix and a vector
def singleRotation(matrix, vector):
    res = [0, 0, 0]
    res[0] = vector[0] * matrix[0][0] + vector[1] * matrix[0][1]\
     + vector[2] * matrix[0][2]
    res[1] = vector[0] * matrix[1][0] + vector[1] * matrix[1][1]\
     + vector[2] * matrix[1][2]
    res[2] = vector[0] * matrix[2][0] + vector[1] * matrix[2][1]\
     + vector[2] * matrix[2][2]
    return res

def rotateVectorX(vector, rx):
    matrixx = [[1, 0, 0],
               [0, cos(rx), sin(rx)],
               [0, -sin(rx), cos(rx)]]
    return singleRotation(matrixx, vector)

def rotateVectorY(vector, ry):
    matrixy = [[cos(ry), 0, -sin(ry)],
               [0, 1, 0],
               [sin(ry), 0, cos(ry)]]
    return singleRotation(matrixy, vector)

def rotateVectorZ(vector, rz):
    matrixz = [[cos(rz), sin(rz), 0],
               [-sin(rz), cos(rz), 0],
               [0, 0, 1]]
    return singleRotation(matrixz, vector)

# Rotates a vector by ry degrees around the y-axis, rx degrees around the x-axis
# and rz degrees around the z-axis. Negative angles mean clockwise
def rotateVector(vector, ry, rx, rz):
    matrixy = [[cos(ry), 0, -sin(ry)],
               [0, 1, 0],
               [sin(ry), 0, cos(ry)]]
    rotated = singleRotation(matrixy, vector)
    matrixx = [[1, 0, 0],
               [0, cos(rx), sin(rx)],
               [0, -sin(rx), cos(rx)]]
    rotated = singleRotation(matrixx, rotated)
    matrixz = [[cos(rz), sin(rz), 0],
               [-sin(rz), cos(rz), 0],
               [0, 0, 1]]
    rotated = singleRotation(matrixz, rotated)

    return rotated

k = 0.01720209894

# Converts a time in solar days to Gaussian days given a reference point
def timeToGaussian(t, t0=0):
    return k * (t - t0)
