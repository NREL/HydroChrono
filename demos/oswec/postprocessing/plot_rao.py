import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re  # import regular expressions module
from scipy.fft import fft

results_directory = '..\..\..\..\HydroChrono_build\demos\RelWithDebInfo\\results\oswec\\regular_waves'

filenames = [fn for fn in os.listdir(results_directory) if fn.startswith('oswec_reg_waves_') and fn.endswith('.txt') and not fn.endswith('_duration.txt')]

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

filenames.sort(key=natural_keys)

num_plots = len(filenames)  # Adjust the number of subplots according to the number of files
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
        axs[i].plot(df['Time (s)'], df['Pitch (rads)'])
        axs[i].set_title(f'Wave Number {i+1}')
    else:
        print(f"Skipping file {filename} as it does not contain 'Pitch (rads)' column.")

# Add x/y labels to the last subplot
axs[-1].set_xlabel('Time (s)')
axs[int(len(axs)/2)].set_ylabel('Pitch (rads)')

plt.tight_layout()
plt.show()

# Initiate empty lists for storing wave periods and corresponding response amplitudes
wave_periods = []
response_amplitudes = []

for i, filename in enumerate(filenames):
    try:
        with open(os.path.join(results_directory, filename), 'r') as file:
            lines = file.readlines()
            # Read the wave amplitude and frequency and convert to wave period in seconds
            wave_amplitude = float(lines[1].split('\t')[1])
            wave_omega = float(lines[2].split('\t')[1])
            print(wave_omega)
            wave_period = 2 * np.pi / wave_omega  # omega = 2pi/T
            wave_periods.append(wave_period)

        # Read time series data
        df = pd.read_csv(os.path.join(results_directory, filename), skiprows=4, delim_whitespace=True, header=None, names=['Time (s)', 'Pitch (rads)'], dtype={'Time (s)': float, 'Pitch (rads)': float})
    except pd.errors.EmptyDataError:
        print(f"Skipping file {filename} as it's empty or unreadable.")
        continue

    if 'Pitch (rads)' in df.columns:
        # Calculate amplitude of the heave response
        length = len(df['Pitch (rads)'])
        response_amplitude = np.abs(np.max(df['Pitch (rads)'][int(0.8*length):]))  # the resting position is -2
        # Calculate RAO by dividing response amplitude by the wave amplitude
        RAO = response_amplitude / wave_amplitude
        response_amplitudes.append(RAO)
    else:
        print(f"Skipping file {filename} as it does not contain 'Pitch (rads)' column.")

# Plotting the RAO
plt.figure(figsize=(10,6))

# Plotting the RAO comparison data
# plt.plot(rao_data['Wave Period (s)'], rao_data['InWave-HOTINT'], 'x--', label='InWave-HOTINT')
# plt.plot(rao_data['Wave Period (s)'], rao_data['InWave (Lin)'], 'x--', label='InWave (Lin)')
# plt.plot(rao_data['Wave Period (s)'], rao_data['Marin (NLin)'], 'x--', label='Marin (NLin)')
# plt.plot(rao_data['Wave Period (s)'], rao_data['NREL (NLin)'], 'x--', label='NREL (NLin)')

plt.plot(wave_periods, response_amplitudes, 'o-', label='Simulation Data', color='black')
plt.plot(rao_data['Wave Period (s)'], rao_data['RAO']*2, 'x--', label='WEC-Sim')

plt.xlabel('Wave Period (s)')
plt.ylabel('RAO (rads/m)')
plt.title('Response Amplitude Operator (RAO)')
plt.legend()
plt.grid(True)
plt.show()
