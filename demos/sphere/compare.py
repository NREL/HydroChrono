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



    testData[:,1] += 2.0  # shift to recenter decay at sea level 0

    nval = testData.shape[0]

    # resample refData to testData sampling rate
    # Suppose dt is constant with same computation physical time  
    x = np.linspace(testData[0, 0], testData[nval-1, 0], nval)
    y1 = np.interp(x, refData[:,0], refData[:,6])
    y2 = np.interp(x, testData[:,0], testData[:,1])

    yd = y1 - y2

    n1 = np.linalg.norm(yd)/nval
    n2 = np.linalg.norm(yd, np.inf)
    
    # Very specific to sphere_decay
    if (n1 > 1e-4 or n2 > 0.02): sys.exit(1)  # error


    #plt.plot(x, y1, '+')
    #plt.plot(x, y2, '--')    
    #plt.plot(x, yd, '--') 
    #plt.show()


    sys.exit(0)