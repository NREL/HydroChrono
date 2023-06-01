import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

# Set some parameters
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# Read HydroChrono results data
decayTestFile = '../../../../HydroChrono_build/demos/RelWithDebInfo/results/sphere_decay.txt'
decayTestData = np.loadtxt(decayTestFile, skiprows=1)

# Read reference data
refFile = './sphere_decay_comparison.xlsx'
refData = pd.read_excel(refFile, 'sphere_decay_comparison')

# Create figure & plot the data
fig, ax = plt.subplots(figsize=(10, 6))  # Set the figure size for a larger, more visible plot
ax.plot(refData['Time (s)'], refData['NREL (CFD)'], '--', label='NREL (CFD)', color='royalblue')
ax.plot(refData['Time (s)'], refData['WavEC (Lin)'], '--', label='WavEC (Lin)', color='forestgreen')
ax.plot(refData['Time (s)'], refData['InWave+H (Lin)'], '--', label='InWave-HOTINT (Lin)', color='darkorange')
ax.plot(decayTestData[:,0], decayTestData[:,1]+2, '--', label='HydroChrono', color='black', linewidth=1.5)

# Set x and y axis limits
ax.set_xlim([0, 40])
ax.set_ylim([-1.0, 1.0])

# Set labels, title & font
ax.set_xlabel('Time (s)', fontsize=14, fontname='Times New Roman')
ax.set_ylabel('Heave position (m)', fontsize=14, fontname='Times New Roman')
ax.set_title('Decay test comparison', fontsize=18, fontname='Times New Roman')

# Set grid, legend, and tight layout
ax.grid(True)
ax.legend(fontsize=12)
plt.tight_layout()

# Save the figure
plt.savefig('sphere_decay_1m.png', dpi=300)
plt.show()