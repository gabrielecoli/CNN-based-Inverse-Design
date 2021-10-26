import numpy as np
import os
import sys

densities = [round(a,2) for a in np.linspace(0.75,1.05,31)]
temperatures = [round(a,2) for a in np.linspace(0.1,0.25,16)]
trials = [0]

for d in densities:
	os.system("mkdir density" + str(d))
	for t in temperatures:
		os.system("mkdir density" + str(d) + "/temperature" + str(t))
		for i in trials:
			os.system("mkdir density" + str(d) + "/temperature" + str(t) + "/trial" + str(i))
			os.system("cp -r StartingPoint/* density" + str(d) + "/temperature" + str(t) + "/trial" + str(i) + "/.")
			fileName = "density" + str(d) + "/temperature" + str(t) + "/trial" + str(i) + "/siminfo.dat"
			fileIn = open(fileName,'w')
			fileIn.write(str(t) + ' ' + str(d) + ' 1.35\n')

print('Done!')
