isImported = True
try: 
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas
    import h5py
except ImportError:
    print("Can't genereate plots because of missing python modules")
    isImported = False
    pass

if isImported:

    #df1 = pandas.read_csv('user/_validation/demos/results/sphere_decay_comparison.txt', sep='\s+', header=None)
    df1 = pandas.read_csv('results/sphere_decay_comparison.txt', sep='\s+', header=None)

    # assert df1.ndim == 2
    # assert df1.shape[1] == 2

    data1 = df1.to_numpy()


    fig = plt.figure(figsize=(8, 6), dpi=500)

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.plot(data1[:,0], data1[:,1], 'ko', markevery=16, markersize=3.5, label='Analytical')
    #ax.plot(data2[:,4], data2[:,5], 'g', label='LiBat3D')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Decay (m)')
    ax.legend()
    ax.grid()
    plt.title(r'Decay')
    plt.show()

