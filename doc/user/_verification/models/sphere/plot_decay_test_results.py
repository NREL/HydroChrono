isImported = True
try: 
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import h5py
except ImportError:
    print("Can't genereate plots because of missing python modules")
    isImported = False
    pass

if isImported:
    # Read the data from the TXT file with multiple spaces as delimiter
    data_dir = '../../../../../demos/sphere/postprocessing/'
    hc_results = pd.read_csv(f'{data_dir}sphere_decay_hc_data.txt', delim_whitespace=True)
    ref_data = pd.read_csv(f'{data_dir}sphere_decay_ref_data.txt', delim_whitespace=True)  # Assuming the file uses tab separation

    # 16:9 aspect ratio
    plt.figure(figsize=(16*0.75, 9*0.75), dpi=500)

    # Read and convert data
    data_hc = hc_results.to_numpy()
    data_ref = ref_data.to_numpy()

    # Plot the data
    plt.plot(data_hc[:, 0], data_hc[:, 1] + 2.0, 'k-', linewidth=2.5, label='HydroChrono')

    # Plot additional data series
    plt.plot(data_ref[:, 0], data_ref[:, 1], '--', linewidth=1.5, label='ProteusDS (Lin)')
    plt.plot(data_ref[:, 0], data_ref[:, 3], '--', linewidth=1.5, label='Marin (Lin)')
    plt.plot(data_ref[:, 0], data_ref[:, 4], '--', linewidth=1.5, label='NREL (CFD)')

    # Labeling and aesthetics
    plt.xlabel(r'Time ($s$)', fontsize=16)
    plt.ylabel(r'Decay ($m$)', fontsize=16)
    plt.title(r'Sphere Decay Test Results', fontsize=18)
    # Customize tick marks and grid lines
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend(fontsize=14)

    plt.tight_layout()

    # Save the plot as a high-resolution PNG image
    plt.savefig("sphere_decay_1m_verification.png", dpi=500, format='png')