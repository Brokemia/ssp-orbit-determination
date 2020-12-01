import odlib
from math import sqrt, acos, sin, cos

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

# Calculates and displays the orbital elements of an asteroid given
# its position at 3 different times
def odMain():
    perihelionTime = 2456968.3929589586
    pos1 = [0.3856417122137088, -1.227683289526401, 0.4709645930434538]
    t1 = 2458312.5
    tau1 = odlib.timeToGaussian(t1, perihelionTime)
    pos2 = [0.3856814471358135, -1.227674472241950, 0.4709778048877576]
    t2 = 2458312.503472222
    tau2 = odlib.timeToGaussian(t2, perihelionTime)
    pos3 = [0.3857211815231323, -1.227665653255291, 0.4709910160789988]
    t3 = 2458312.506944444
    tau3 = odlib.timeToGaussian(t3, perihelionTime)

    drx = (pos3[0] - pos1[0]) / (tau3 - tau1)
    dry = (pos3[1] - pos1[1]) / (tau3 - tau1)
    drz = (pos3[2] - pos1[2]) / (tau3 - tau1)

    rDot = [drx, dry, drz]

    elements = orbitalElements(pos2, rDot)
    realElements = [1.056800391773215, .3442329516328222, 25.15526140011037,
            236.2379969538250, 255.5046626217337, 139.5152750732308]

    elements[2] = odlib.radiansToDegrees(elements[2])
    elements[3] = odlib.radiansToDegrees(elements[3])
    elements[4] = odlib.radiansToDegrees(elements[4])
    elements[5] = odlib.radiansToDegrees(elements[5])

    while elements[4] < 0:
        elements[4] += 360
    while elements[4] > 360:
        elements[4] -= 360

    percentErrors = []
    for i in range(len(elements)):
        # I know this is over 80 characters
        # but I can't really make it shorter without shortening variable names
        percentErrors.append(100 * abs(elements[i] - realElements[i]) / realElements[i])

    print("Name\tExpected Value\t\tCalculated Value\tPercent Error")
    print(f"a\t{realElements[0]}\t{elements[0]}\t{percentErrors[0]}%")
    print(f"e\t{realElements[1]}\t{elements[1]}\t{percentErrors[1]}%")
    print(f"i\t{realElements[2]}\t{elements[2]}\t{percentErrors[2]}%")
    print(f"Omega\t{realElements[3]}\t{elements[3]}\t{percentErrors[3]}%")
    print(f"w\t{realElements[4]}\t{elements[4]}\t{percentErrors[4]}%")
    print(f"M\t{realElements[5]}\t{elements[5]}\t{percentErrors[5]}%")

    return elements

odMain()
