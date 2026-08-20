import ezodf
import numpy as np


def handle_spectrumnpz_data(filename, frequency_key="arr_0", amplitude_key="arr_1", **kwargs):
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
    spreadsheet = ezodf.newdoc(doctype="ods", filename=file_name)
    max_len = max(len(arr) for arr in arrays)
    sheet = ezodf.Sheet(sheet_name, size=(max_len + 1, len(arrays)))
    spreadsheet.sheets += sheet

    for i, arr in enumerate(arrays):
        sheet[(0, i)].set_value(str(names[i]))
        for j in range(len(arr)):
            sheet[j + 1, i].set_value(arr[j])

    # Save the spreadsheet
    spreadsheet.save()
