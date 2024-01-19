import numpy as np
import sys

# Function to read two wave data files
def read_wave_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    wave_data = {}

    heave_data = []
    for line in lines[1:]:
        time, heave = line.split()
        heave_data.append(float(heave))
    wave_data['heave'] = np.array(heave_data)

    return wave_data

def calculate_l2_norm(heave_data1, heave_data2):
    yd = heave_data1 - heave_data2
    nval = yd.shape[0]
    l2 = np.linalg.norm(yd)/nval
    return l2

def compare_wave_files(file_path1, file_path2):
    wave_data1 = read_wave_data(file_path1)
    wave_data2 = read_wave_data(file_path2)

    # calculate the l2 norm of the difference in heave
    diff_heave = calculate_l2_norm(wave_data1['heave'],wave_data2['heave'])

    # If the difference is more than 1e-6 then error out
    if (diff_heave > 1e-6): sys.exit(1)  # error



if __name__ == "__main__":

    fname_ref = sys.argv[1]
    fname_rst = sys.argv[2]

    compare_wave_files(fname_ref, fname_rst)
