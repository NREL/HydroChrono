import sys
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    """
    Compare result from reference file
    
    Usage: >
    """

    fname_ref = sys.argv[1]
    fname_rst = sys.argv[2]

    print(fname_ref)

    refData = np.loadtxt(fname_ref, skiprows=1)

    testData = np.loadtxt(fname_rst, skiprows=1)

    print(refData.shape)
    print(testData.shape)


    # Don't think this is needed for this demo
    # testData[:,1] += 2.0  # shift to recenter decay at sea level 0

    nval = testData.shape[0]

    # resample refData to testData sampling rate
    # Suppose dt is constant with same computation physical time  
    x = np.linspace(testData[0, 0], testData[nval-1, 0], nval)
    # Comapre the flap's pitch in radian
    flapPitchRef = np.interp(x, refData[:,0], refData[:,1])
    flapPitchTest = np.interp(x, testData[:,0], testData[:,1])
    flapPitchComp = flapPitchRef - flapPitchTest


    #Forbnius norm - Flap pitch
    flapPitchn1 = np.linalg.norm(flapPitchComp)/nval
    #infinity norm - Flap pitch
    flapPitchn2 = np.linalg.norm(flapPitchComp, np.inf)

    
    # Needs to be changed based on domain knowledge, for now this is random
    if (flapPitchn1 > 1e-4 or flapPitchn2 > 0.02): sys.exit(1)  # error


    #plt.plot(x, y1, '+')
    #plt.plot(x, y2, '--')    
    #plt.plot(x, yd, '--') 
    #plt.show()


    sys.exit(0)