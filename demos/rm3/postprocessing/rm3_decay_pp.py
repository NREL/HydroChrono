import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

decayDataPath = f'..\\..\\..\\..\\HydroChrono_build\\Release\\results\\rm3\\decay\\rm3_decay.txt'
decayData = np.loadtxt(decayDataPath, skiprows=1)

refFile = './rm3_ref_data.xlsx'
refData = pd.read_excel(refFile, '0.01mDecay_0ptoDamping')

fig, axs = plt.subplots(2, 1)
axs[0].plot(refData['Time (s)'], refData['Float Heave (m)'], label='wec-sim')
axs[0].plot(decayData[:,0], decayData[:,1], '--', label='chrono')
axs[0].set_xlim([0.0, 30.0])
axs[0].set_xlabel('Time (s)')
axs[0].legend()

axs[1].plot(refData['Time (s)'], refData['Spar Heave (m)'], label='wec-sim')
axs[1].plot(decayData[:,0], decayData[:,2], '--', label='chrono')
axs[1].set_xlim([0.0, 30.0])
axs[1].set_xlabel('Time (s)')
axs[1].legend()

fig.tight_layout()
plt.show()