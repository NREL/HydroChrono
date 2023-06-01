import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import os
import re
from scipy.fft import fft

# Set some parameters
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
results_directory = '..\..\..\..\HydroChrono_build\demos\RelWithDebInfo\\results\oswec\\regular_waves'

filenames = [fn for fn in os.listdir(results_directory) if fn.startswith('oswec_reg_waves_') and fn.endswith('.txt') and not fn.endswith('_duration.txt')]

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

filenames.sort(key=natural_keys)

num_plots = len(filenames)
fig, axs = plt.subplots(num_plots, figsize=(10, 2 * num_plots))

# Read the RAO comparison data
rao_data = pd.read_excel('RAO_dat.xlsx', sheet_name='PTO_damping=0.0')

for i, filename in enumerate(filenames):
    try:
        df = pd.read_csv(os.path.join(results_directory, filename), skiprows=4, delim_whitespace=True, header=None, names=['Time (s)', 'Pitch (rads)'], dtype={'Time (s)': float, 'Pitch (rads)': float})
    except pd.errors.EmptyDataError:
        print(f"Skipping file {filename} as it's empty or unreadable.")
        continue

    if 'Pitch (rads)' in df.columns:
        axs[i].plot(df['Time (s)'], df['Pitch (rads)'], '-', color='black', linewidth=1.5)
        axs[i].set_title(f'Wave Number {i+1}', fontsize=18, fontname='Times New Roman')
    else:
        print(f"Skipping file {filename} as it does not contain 'Pitch (rads)' column.")

# Add x/y labels to the last subplot
axs[-1].set_xlabel('Time (s)', fontsize=14, fontname='Times New Roman')
axs[int(len(axs)/2)].set_ylabel('Pitch (rads)', fontsize=14, fontname='Times New Roman')

for ax in axs:
    ax.grid(True)

plt.tight_layout()
plt.show()

wave_periods = []
response_amplitudes = []

for i, filename in enumerate(filenames):
    try:
        with open(os.path.join(results_directory, filename), 'r') as file:
            lines = file.readlines()
            wave_amplitude = float(lines[1].split('\t')[1])
            wave_omega = float(lines[2].split('\t')[1])
            wave_period = 2 * np.pi / wave_omega  # omega = 2pi/T
            wave_periods.append(wave_period)

        df = pd.read_csv(os.path.join(results_directory, filename), skiprows=4, delim_whitespace=True, header=None, names=['Time (s)', 'Pitch (rads)'], dtype={'Time (s)': float, 'Pitch (rads)': float})
    except pd.errors.EmptyDataError:
        print(f"Skipping file {filename} as it's empty or unreadable.")
        continue

    if 'Pitch (rads)' in df.columns:
        length = len(df['Pitch (rads)'])
        response_amplitude = np.abs(np.max(df['Pitch (rads)'][int(0.8*length):]))
        RAO = response_amplitude / wave_amplitude
        response_amplitudes.append(RAO)
    else:
        print(f"Skipping file {filename} as it does not contain 'Pitch (rads)' column.")

# Plotting the RAO
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(wave_periods, response_amplitudes, 'o-', label='HydroChrono', color='black', linewidth=1.5)
ax.plot(rao_data['Wave Period (s)'], rao_data['RAO']*2, 'x--', label='WEC-Sim', color='royalblue')

ax.set_xlabel('Wave Period (s)', fontsize=14, fontname='Times New Roman')
ax.set_ylabel('RAO (rads/m)', fontsize=14, fontname='Times New Roman')
ax.set_title('Response Amplitude Operator (RAO)', fontsize=18, fontname='Times New Roman')

ax.grid(True)
ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig('RAO_comparison.png', dpi=300)
plt.show()
