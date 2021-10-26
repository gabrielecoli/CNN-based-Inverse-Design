import numpy as np
import os
import sys

densities = [round(a,2) for a in np.linspace(0.75,1.05,31)]
temperatures = [round(a,2) for a in np.linspace(0.1,0.25,16)]
trials = [0]

for d in densities:
        for t in temperatures:
                for i in trials:
                        os.chdir("density" + str(d) + "/temperature" + str(t) + "/trial" + str(i))
                        os.system("qsub -V -cwd -q all.q@node[0-5][0-9].cm.cluster run_cnn.sh")
			os.chdir("../../../")

print('Done!')
