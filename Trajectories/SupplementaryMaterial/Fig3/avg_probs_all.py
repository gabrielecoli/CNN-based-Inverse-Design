import numpy as np
import os
import sys

densities = [round(a,2) for a in np.linspace(0.75,1.05,31)]
temperatures = [round(a,2) for a in np.linspace(0.1,0.25,16)]
trials = [0]
PHASES = ['fluid','square','hexagonal','qc12','qc10','qc18']
N_PHASES = len(PHASES)
N_TRIALS = len(trials)
cols = [i for i in range(N_PHASES)]

fileOut = []
for phase in PHASES:
	fileNameOut = 'AVG/CNN/' + phase + '.dat'
	fileOut.append(open(fileNameOut,'w'))

for d in densities:
	print(d)
	for t in temperatures:
		shape = (N_PHASES,N_TRIALS)
		all_probs = np.zeros(shape)
		for i in trials:
			filePath = "density" + str(d) + "/temperature" + str(t) + "/trial" + str(i) + "/"
			fileNameInfo = filePath + "siminfo.dat"
			fileNameCNN = filePath + "probs_cnn.dat"
			if(i==trials[0]):
				thisInfo = np.loadtxt(fileNameInfo)
				print(thisInfo)
			all_probs[:,0] = np.loadtxt(fileNameCNN,usecols=cols)
		avg_probs = np.mean(all_probs,axis=1)
		for i,phase in enumerate(PHASES):
			fileOut[i].write(str(thisInfo[0]) + ' ' + str(thisInfo[1]) + ' ' + str(round(avg_probs[i],5)) + '\n')


for i in range(N_PHASES):
        fileOut[i].close()



print('Done!')
