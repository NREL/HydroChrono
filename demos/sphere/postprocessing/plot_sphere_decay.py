import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

decayTestFile = '../../../../HydroChrono_build/Release/sphere_decay.txt'
decayTestData = np.loadtxt(decayTestFile, skiprows=1)

refFile = './sphere_decay_comparison.xlsx'
refData = pd.read_excel(refFile, '1m (2)')

plt.plot(refData['Time (s)'], refData['NREL (CFD)'], '-')
plt.plot(decayTestData[:,0], decayTestData[:,1]+2, '--')
plt.show()