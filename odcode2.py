# The only reason this program takes so long is because the terminal
# prints slowly
xbus = 0
xbm = 0
vbus = 65
vbm = 85
time = 0
timeStep = .00001
tiny = .00001
iterations = 0
maxiteration = 1000000

print("Iteration\tt\t\t\tXbm\tXbus")
# time < 6/60 check included to stop the batmobile from catching up while
# on the starting line
while (iterations < maxiteration and abs(xbus - xbm) > tiny) or time < 6/60:
    xbus = vbus * time
    xbm = max(0, vbm * (time - 6/60))

    time += timeStep
    iterations += 1
    # Time converted to minutes
    print(f"{iterations}\t\t{time*60}\t{xbm}\t{xbus}")

if iterations >= maxiteration:
    print("Error: Solution not found within maximum iterations")
else:
    print("The batmobile caught up in", time*60, "minutes.")
