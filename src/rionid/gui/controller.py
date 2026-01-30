from numpy import argsort, where, append
from loguru import logger
from rionid.core import ImportData
from rionid.external.barion.amedata import AMEData

def import_controller(datafile=None, filep=None, alphap=None, refion=None, harmonics=None, 
                      nions=None, amplitude=None, circumference=None, mode=None, value=None, 
                      reload_data=None, remove_baseline=False, psd_baseline_removed_l=1e6,
                      peak_threshold_pct=0.05, min_distance=10, highlight_ions=None, 
                      io_params=None, sim_scalingfactor=None, matching_freq_min=None, 
                      matching_freq_max=None, correct=None):
    """
    Main orchestration function for the RionID simulation workflow.

    This function acts as the bridge between the GUI/CLI inputs and the physics core.
    It initializes the `ImportData` model, loads experimental data, performs the 
    physics calculations (mass-to-charge, revolution frequencies), runs the 
    simulation, and saves the results to disk.

    Parameters
    ----------
    datafile : str, optional
        Path to the experimental data file (.npz, .csv, .root, etc.).
    filep : str, optional
        Path to the particle list file (e.g., LISE++ .lpp output).
    alphap : float, optional
        Momentum compaction factor of the ring. If the value provided is > 1, 
        it is treated as Gamma Transition (gamma_t) and converted automatically 
        via `alphap = 1 / gamma_t^2`.
    refion : str, optional
        The reference ion string in the format 'AAEl+QQ' (e.g., '72Ge+35').
    harmonics : str or list of float, optional
        Harmonic numbers to simulate. Can be a list or a space/comma-separated string.
    nions : int, optional
        If provided, filters the output to show only the top `nions` species 
        sorted by yield (plus the reference ion).
    amplitude : int, optional
        Display option for spectral lines (0 for constant height, 1 for scaled).
    circumference : float, optional
        Ring circumference in meters.
    mode : str, optional
        The operation mode. Options: 'Frequency', 'Brho', 'Kinetic Energy', 'Gamma'.
    value : float, optional
        The numerical value corresponding to the selected `mode` (e.g., the 
        reference frequency in Hz if mode is 'Frequency').
    reload_data : bool, optional
        If True, reloads experimental data from the raw file; otherwise loads 
        from the cached .npz.
    remove_baseline : bool, optional
        If True, applies baseline subtraction algorithms to the experimental spectrum.
    psd_baseline_removed_l : float, optional
        The smoothness parameter (lambda) for the baseline subtraction algorithm (BrPLS).
        Default is 1e6.
    peak_threshold_pct : float, optional
        Relative threshold for peak detection (0.0 to 1.0). Default is 0.05 (5%).
    min_distance : float, optional
        Minimum distance (in data points) between detected peaks.
    highlight_ions : str, optional
        Comma-separated string of ion names to highlight in the plot (e.g., '72Ge+35, 74Se+34').
    io_params : dict, optional
        Dictionary containing specific I/O parameters (e.g., keys for NPZ files, 
        histogram names for ROOT files).
    sim_scalingfactor : float, optional
        Factor to scale the simulated yield/amplitude.
    matching_freq_min : float, optional
        Minimum frequency (Hz) bound for the peak matching algorithm.
    matching_freq_max : float, optional
        Maximum frequency (Hz) bound for the peak matching algorithm.
    correct : list of float, optional
        Coefficients [a0, a1, a2] for a second-order polynomial correction 
        applied to the simulated frequencies.

    Returns
    -------
    ImportData
        The populated data model object containing simulation results and 
        experimental data.
    Exception
        Returns the exception object if an error occurs during processing.
    """
    try:
        # initializations
        if float(alphap) > 1: alphap = 1/float(alphap)**2 # handling alphap and gammat
        fref = brho = ke = gam = None
        if mode == 'Frequency': fref = float(value)
        elif mode == 'Brho': brho = float(value)
        elif mode == 'Kinetic Energy': ke = float(value)
        elif mode == 'Gamma': gam = float(value)
        # Calculations 
        mydata = ImportData(refion, float(alphap), filename=datafile, reload_data=reload_data, 
                            circumference=circumference, highlight_ions=highlight_ions,
                            remove_baseline=remove_baseline, psd_baseline_removed_l=psd_baseline_removed_l,
                            peak_threshold_pct=peak_threshold_pct, min_distance=min_distance,
                            matching_freq_min=matching_freq_min, matching_freq_max=matching_freq_max,
                            io_params=io_params)
        
        mydata._set_particles_to_simulate_from_file(filep)
        mydata._calculate_moqs()
        mydata._calculate_srrf(fref=fref, brho=brho, ke=ke, gam=gam, correct=correct)
        
        if isinstance(harmonics, str):
            harmonics = [float(h.strip()) for h in harmonics.replace(',', ' ').split()]
        elif not isinstance(harmonics, list):
            harmonics = [float(harmonics)]
            
        mydata._simulated_data(brho=brho, harmonics=harmonics, mode=mode, 
                               sim_scalingfactor=sim_scalingfactor, nions=nions)
        # "Outputs"
        if nions:
            display_nions(int(nions), mydata.yield_data, mydata.nuclei_names, mydata.simulated_data_dict, refion, harmonics)
        
        logger.info(f'Simulation results (ordered by frequency) will be saved to simulation_result.out')
        sort_index = argsort(mydata.srrf)
        # Save the results to a file with the specified format
        save_simulation_results(mydata, mode, harmonics, sort_index)
        logger.info(f'Succesfully saved!')

        return mydata # Returns the simulated spectrum data 
    
    except Exception as e:
        print(f"Error during calculations: {str(e)}")
        raise e

def display_nions(nions, yield_data, nuclei_names, simulated_data_dict, ref_ion, harmonics):
    """
    Filters the simulated data to retain only the top N ions based on yield.

    This function modifies the `simulated_data_dict` in-place. It ensures that
    the Reference Ion is always included in the list, even if its yield is low.

    Parameters
    ----------
    nions : int
        The number of ions to keep.
    yield_data : array-like
        Array of yield values for all simulated ions.
    nuclei_names : array-like
        Array of ion names corresponding to the yield data.
    simulated_data_dict : dict
        Dictionary mapping harmonic numbers to simulation arrays.
    ref_ion : str
        The name of the reference ion (e.g., '72Ge+35').
    harmonics : list of float
        The list of harmonics being simulated.
    """
    # 1. Get indices of the top N ions
    sorted_indices = argsort(yield_data)[::-1][:nions]
    
    # 2. Find the index of the reference ion
    ref_index = where(nuclei_names == ref_ion)[0]
    
    # 3. FIX: Check if ref_index is not empty AND if the scalar index is missing
    if ref_index.size > 0 and ref_index[0] not in sorted_indices:
        sorted_indices = append(sorted_indices, ref_index[0])
        
    # 4. Filter the names
    nuclei_names = nuclei_names[sorted_indices]
    
    # 5. Filter the simulated data for each harmonic
    for harmonic in harmonics: 
        name = f'{harmonic}'
        if name in simulated_data_dict:
            simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]

def save_simulation_results(mydata, mode, harmonics, sort_index, filename='simulation_result.out'):
    """
    Saves the simulation results to a formatted text file.

    Writes a table containing Ion Name, Frequency, Yield, m/q, and Mass for
    every simulated ion, sorted according to `sort_index`.

    Parameters
    ----------
    mydata : ImportData
        The data object containing the calculated physics results (srrf, moq, etc.).
    mode : str
        The simulation mode ('Frequency', 'Brho', etc.) used to determine how
        absolute frequencies are calculated.
    harmonics : list of float
        The list of harmonics used in the simulation.
    sort_index : array-like
        Indices used to sort the output (typically sorted by frequency).
    filename : str, optional
        The output filename. Default is 'simulation_result.out'.
    """
    with open(filename, 'w') as file:
        # Writing harmonics and brho information
        brho = mydata.brho
        for harmonic in harmonics:
            header0 = f'Harmonic: {harmonic} , Bp: {brho:.6f} [Tm]'
            logger.info(header0)
            file.write(header0 + '\n')
        
        # Writing the header for the data table
        header1 = f"{'ion':<15}{'fre[Hz]':<30}{'yield [pps]':<15}{'m/q [u]':<15}{'m [eV]':<15}"
        file.write(header1 + '\n')
        file.write('-' * len(header1) + '\n')
        logger.info(header1)
        
        # Writing the sorted simulation results
        for i in sort_index:
            ion = mydata.nuclei_names[i]
            if mode == 'Frequency': fre = mydata.srrf[i] * mydata.ref_frequency
            elif mode == 'Brho': fre = mydata.srrf[i] * mydata.ref_frequency*harmonic
            yield_ = mydata.yield_data[i]
            moq = mydata.moq[ion]
            mass_u = mydata.total_mass[ion]
            mass = AMEData.to_mev(mass_u) * 1e6
            result_line = f"{ion:<15}{fre:<30.10f}{yield_:<15.4e}{moq:<15.12f}{mass:<15.3f}"
            logger.info(result_line)
            file.write(result_line + '\n')