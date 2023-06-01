import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

# Set some parameters
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

# Read HydroChrono results data
decayTestFile = '../../../../HydroChrono_build/demos/RelWithDebInfo/results/oswec_decay.txt'
decayTestData = np.loadtxt(decayTestFile, skiprows=1)

# Read reference data
refFile = './wecsim_oswec_decay.txt'
refData = np.loadtxt(refFile, skiprows=1)

# Create figure & plot the data
fig, ax = plt.subplots(figsize=(10, 6))  # Set the figure size for a larger, more visible plot
ax.plot(decayTestData[:,0], decayTestData[:,1], '-', label='HydroChrono', color='black', linewidth=1.5)
ax.plot(refData[:,0], refData[:,1], '--', label='WEC-Sim', color='royalblue')

# Set x and y axis limits
ax.set_xlim([0, 400])
# ax.set_ylim([-1.0, 1.0])

# Set labels, title & font
ax.set_xlabel('Time (s)', fontsize=14, fontname='Times New Roman')
ax.set_ylabel('Pitch angle (rads)', fontsize=14, fontname='Times New Roman')
ax.set_title('Decay test comparison', fontsize=18, fontname='Times New Roman')

# Set grid, legend, and tight layout
ax.grid(True)
ax.legend(fontsize=12)
plt.tight_layout()

# Save the figure
plt.savefig('oswec_decay_pitch.png', dpi=300)
plt.show()