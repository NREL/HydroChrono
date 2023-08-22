# Compare to Chrono itself from commit 745bdabd60c

import numpy as np
import sys
# Function to read two wave data files
def read_wave_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)
    return data

def calculate_l2_norm(data1, data2):
    yd = data1 - data2
    nval = yd.shape[0]
    l2 = np.linalg.norm(yd)/nval
    return l2

def compare_wave_files(file_path1, file_path2):
    wave_data1 = read_wave_data(file_path1)
    wave_data2 = read_wave_data(file_path2)

    # calculate the l2 norm of the difference in the base surge
    diff_base_surge = calculate_l2_norm(wave_data1[:,1],wave_data2[:,1])
    # If the difference is more than 1e-6 then error out
    if (diff_base_surge > 1e-6): sys.exit(1)  # error

    # calculate the difference in the base pitch
    diff_base_pitch = calculate_l2_norm(wave_data1[:,2], wave_data2[:,2])
    # If the difference is more than 1e-10 then error out (since radian scale is smaller a smaller value is used here)
    if (diff_base_pitch > 1e-10): sys.exit(1)  # error
    
    # calculate the difference in the flap fore pitch 
    diff_flap_fore_pitch = calculate_l2_norm(wave_data1[:,3], wave_data2[:,3])

    if (diff_flap_fore_pitch > 1e-6): sys.exit(1)  # error

    # calculate the difference in the flap aft pitch 
    diff_flap_aft_pitch = calculate_l2_norm(wave_data1[:,4], wave_data2[:,4])

    if (diff_flap_aft_pitch > 1e-6): sys.exit(1)  # error


if __name__ == "__main__":

    fname_ref = sys.argv[1]
    fname_rst = sys.argv[2]

    compare_wave_files(fname_ref, fname_rst)