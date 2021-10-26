import numpy as np
import os
import sys

densities = [round(a,2) for a in np.linspace(0.75,1.05,31)]
temperatures = [round(a,2) for a in np.linspace(0.1,0.25,16)]
trials = [0]
N_CHI = 24
N_TRIALS = len(trials)
cols = [i+2 for i in range(N_CHI)]

fileOut = []
for chi in range(N_CHI):
	fileNameOut = 'AVG/CHI/chi' + str(chi+1) + '.dat'
	fileOut.append(open(fileNameOut,'w'))

for d in densities:
	print(d)
	for t in temperatures:
		shape = (N_CHI,N_TRIALS)
		all_chi = np.zeros(shape)
		for i in trials:
			filePath = "density" + str(d) + "/temperature" + str(t) + "/trial" + str(i) + "/"
			fileNameInfo = filePath + "siminfo.dat"
			fileNameChi = filePath + "boop.dat"
			if(i==trials[0]):
				thisInfo = np.loadtxt(fileNameInfo)
				print(thisInfo)
			thisChi = np.loadtxt(fileNameChi,usecols=cols)
			all_chi[:,0] = np.mean(thisChi,axis=0)
		avg_chi = np.mean(all_chi,axis=1)
		for chi in range(N_CHI):
			fileOut[chi].write(str(thisInfo[0]) + ' ' + str(thisInfo[1]) + ' ' + str(round(avg_chi[chi],5)) + '\n')


for chi in range(N_CHI):
        fileOut[chi].close()



print('Done!')
