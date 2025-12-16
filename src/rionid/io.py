from iqtools.tools import read_rsa_specan_xml, read_rsa_data_csv, read_rsa_result_csv
import numpy as np
import os

def read_tdsm_bin(path):
    base_path, _ = os.path.splitext(path)
    bin_fre_path = os.path.join(base_path + '.bin_fre')
    bin_time_path = os.path.join(base_path + '.bin_time')
    bin_amp_path = os.path.join(base_path + '.bin_amp')

    try:
        fre = np.fromfile(bin_fre_path, dtype=np.float64)
        time = np.fromfile(bin_time_path, dtype=np.float32)
        #amp = np.fromfile(bin_amp_path, dtype=np.float32)
        amp = np.memmap(bin_amp_path, dtype=np.float32, mode='r', shape=(len(time), len(fre)))
    except IOError as e:
        raise Exception(f"Error reading files: {e}")
    
    if len(time) == 0 or len(fre) == 0:
        raise ValueError("Time or frequency data files are empty")
    
    try:
        amp = amp.reshape((len(time), len(fre)))
    except ValueError as e:
        raise ValueError(f"Amplitude data cannot be reshaped: {e}")

    midpoint = len(fre) // 2
    frequency = np.concatenate((fre[midpoint:], fre[:midpoint]))
    amplitude = np.concatenate((amp[:, midpoint:], amp[:, :midpoint]), axis=1)

    return frequency, time, amplitude

def handle_read_tdsm_bin(path):
    frequency, _, amplitude = read_tdsm_bin(path)
    amplitude_avg = np.average(amplitude, axis = 0)
    return frequency, amplitude_avg

def handle_read_rsa_specan_xml(filename):
    #frequency, amplitude, _ = read_rsa_specan_xml(filename)
    freq, power, _ = read_rsa_specan_xml(self.filename)
    power = power - power.min() #in order to avoid negative values: power - (-|value_min|) > #power
    #normalized_power = power / power.max()
    return freq, power
    #return frequency, amplitude

def handle_read_rsa_data_csv(filename):
    data = read_rsa_data_csv(filename)
    #obtain frequency, power from data
    return data

def handle_read_rsa_result_csv(filename):
    frequency, amplitude = read_rsa_result_csv(filename)
    return frequency, amplitude

def handle_tiqnpz_data(filename):
    data = np.load(filename)
    frequency = data['arr_0'].flatten()
    amplitude = data['arr_2']
    amplitude_average = np.average(amplitude[5:,:], axis=0)
    return frequency, amplitude_average

def handle_spectrumnpz_data(filename):
    data = np.load(filename)
    frequency = data['arr_0'].flatten()
    amplitude = data['arr_1']
    return frequency, amplitude

def handle_prerionidnpz_data(filename):
    data = np.load(filename)
    frequency = data['x']
    amplitude = data['y']
    return frequency, amplitude

def read_psdata(filename, dbm = False):
    if dbm: 
        frequency, amplitude = np.genfromtxt(filename, skip_header = 1, delimiter='|', usecols = (0,2))
    else: 
        frequency, amplitude = np.genfromtxt(filename, skip_header = 1, delimiter='|', usecols = (0,1))

    return frequency, amplitude

def write_arrays_to_ods(file_name, sheet_name, names, *arrays):
    # Create the ods spreadsheet and add a sheet
    spreadsheet = ezodf.newdoc(doctype='ods', filename=file_name)
    max_len = max(len(arr) for arr in arrays)
    sheet = ezodf.Sheet(sheet_name,size=(max_len+1,len(arrays)))
    spreadsheet.sheets += sheet
    
    for i, arr in enumerate(arrays):
        sheet[(0, i)].set_value(str(names[i]))
        for j in range(len(arr)):
            sheet[j+1, i].set_value(arr[j])

    # Save the spreadsheet
    spreadsheet.save()
