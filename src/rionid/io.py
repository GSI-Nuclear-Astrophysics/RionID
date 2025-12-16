import numpy as np
import ezodf
import os

def read_tdsm_bin(path):
    """
    Reads custom binary TDSM files (.bin_fre, .bin_time, .bin_amp).

    This function expects three files to exist with the same base name but different
    extensions. It uses memory mapping for the amplitude file to handle large datasets efficiently.
    It also performs an FFT shift operation (swapping halves) on the frequency and amplitude arrays.

    Parameters
    ----------
    path : str
        The file path (extension is ignored, base name is used).

    Returns
    -------
    tuple
        (frequency_array, time_array, amplitude_matrix)
    
    Raises
    ------
    Exception
        If files cannot be read or dimensions are inconsistent.
    """
    base_path, _ = os.path.splitext(path)
    bin_fre_path = os.path.join(base_path + '.bin_fre')
    bin_time_path = os.path.join(base_path + '.bin_time')
    bin_amp_path = os.path.join(base_path + '.bin_amp')

    try:
        fre = np.fromfile(bin_fre_path, dtype=np.float64)
        time = np.fromfile(bin_time_path, dtype=np.float32)
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
    """
    Wrapper for reading TDSM binary files and averaging the amplitude.

    Parameters
    ----------
    path : str
        Path to the binary file set.

    Returns
    -------
    tuple
        (frequency_array, averaged_amplitude_array)
    """
    frequency, _, amplitude = read_tdsm_bin(path)
    amplitude_avg = np.average(amplitude, axis=0)
    return frequency, amplitude_avg

def handle_tiqnpz_data(filename, frequency_key='arr_0', amplitude_key='arr_2', time_key='arr_1', **kwargs):
    """
    Handles standard IQTools NPZ files (Time-Frequency-Amplitude).

    This function averages the amplitude over time.
    Note: It skips the first 5 time slices (`amp[5:,:]`) to avoid startup artifacts/shifts.

    Parameters
    ----------
    filename : str
        Path to the .npz file.
    frequency_key : str, optional
        Key for frequency array. Default 'arr_0'.
    amplitude_key : str, optional
        Key for amplitude matrix. Default 'arr_2'.
    time_key : str, optional
        Key for time array. Default 'arr_1'.

    Returns
    -------
    tuple
        (frequency_array, averaged_amplitude_array)
    """
    data = np.load(filename)
    freq = data[frequency_key].flatten()
    amp = data[amplitude_key]
    
    # Averaging from index 5 onwards to remove potential startup artifacts
    amplitude_average = np.average(amp[5:,:], axis=0) 
    return freq, amplitude_average

def handle_spectrumnpz_data(filename, frequency_key='arr_0', amplitude_key='arr_1', **kwargs):
    """
    Handles simple 1D Spectrum NPZ files.

    Parameters
    ----------
    filename : str
        Path to the .npz file.
    frequency_key : str, optional
        Key for frequency array. Default 'arr_0'.
    amplitude_key : str, optional
        Key for amplitude array. Default 'arr_1'.

    Returns
    -------
    tuple
        (frequency_array, amplitude_array)
    """
    data = np.load(filename)
    return data[frequency_key].flatten(), data[amplitude_key]

def handle_prerionidnpz_data(filename):
    """
    Handles legacy PreRionID NPZ files using 'x' and 'y' keys.

    Parameters
    ----------
    filename : str
        Path to the .npz file.

    Returns
    -------
    tuple
        (frequency_array, amplitude_array)
    """
    data = np.load(filename)
    frequency = data['x']
    amplitude = data['y']
    return frequency, amplitude

def read_psdata(filename, dbm=False):
    """
    Reads generic CSV/TXT spectrum files.

    Assumes a pipe-delimited format ('|').

    Parameters
    ----------
    filename : str
        Path to the file.
    dbm : bool, optional
        If True, reads the 3rd column (index 2). If False, reads the 2nd column (index 1).
        Default is False.

    Returns
    -------
    tuple
        (frequency_array, amplitude_array)
    """
    if dbm: 
        frequency, amplitude = np.genfromtxt(filename, skip_header=1, delimiter='|', usecols=(0,2))
    else: 
        frequency, amplitude = np.genfromtxt(filename, skip_header=1, delimiter='|', usecols=(0,1))

    return frequency, amplitude

def write_arrays_to_ods(file_name, sheet_name, names, *arrays):
    """
    Writes data arrays to an OpenDocument Spreadsheet (.ods).

    Parameters
    ----------
    file_name : str
        Output filename.
    sheet_name : str
        Name of the sheet to create.
    names : list of str
        List of column headers.
    *arrays : list of array-like
        Variable number of arrays to write as columns.
    """
    # Create the ods spreadsheet and add a sheet
    spreadsheet = ezodf.newdoc(doctype='ods', filename=file_name)
    max_len = max(len(arr) for arr in arrays)
    sheet = ezodf.Sheet(sheet_name, size=(max_len+1, len(arrays)))
    spreadsheet.sheets += sheet
    
    for i, arr in enumerate(arrays):
        sheet[(0, i)].set_value(str(names[i]))
        for j in range(len(arr)):
            sheet[j+1, i].set_value(arr[j])

    # Save the spreadsheet
    spreadsheet.save()