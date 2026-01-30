import argparse
import os
import logging as log
import sys
from numpy import argsort, where, append, shape
from PyQt5.QtWidgets import QApplication

from rionid.core import ImportData
from rionid.gui.plot import CreatePyGUI
from rionid.io import write_arrays_to_ods

def main():
    """
    Main entry point for the command-line interface (CLI).

    Parses arguments, initializes logging, and dispatches the analysis controller
    for one or more data files.
    """
    scriptname = 'RionID' 
    parser = argparse.ArgumentParser(description="RionID: Ring-stored ion Identification")
    modes = parser.add_mutually_exclusive_group(required=True)

    # Main Arguments
    parser.add_argument('datafile', type=str, nargs='+', help='Name of the input file with data.')
    parser.add_argument('-ap', '--alphap', type=float, help='Momentum compaction factor of the ring.')
    parser.add_argument('-r', '--refion', type=str, help='Reference ion (Format: AAXX+CC, e.g., 72Ge+35).')
    parser.add_argument('-psim', '--filep', type=str, help='Path to particle list file (LISE++ output).')
    parser.add_argument('-hrm', '--harmonics', type=float, default=[1.0], nargs='+', help='Harmonics to simulate.')

    # Secondary Arguments
    parser.add_argument('-n', '--nions', type=int, help='Number of ions to display, sorted by yield.')

    # Arguments for Each Mode (Exclusive)
    modes.add_argument('-b', '--brho', type=float, help='Brho value [Tm] (Isochronous mode).')
    modes.add_argument('-ke', '--kenergy', type=float, help='Kinetic energy [MeV/u] (Isochronous mode).')
    modes.add_argument('-gam', '--gamma', type=float, help='Lorentz factor gamma.')
    modes.add_argument('-f', '--fref', type=float, help='Revolution frequency [Hz] (Standard mode).')
    
    # Visualization & Output
    parser.add_argument('-d', '--ndivs', type=int, default=4, help='Number of divisions (Deprecated).')
    parser.add_argument('-am', '--amplitude', type=int, default=0, help='Display options (0=const, else=scaled).')
    parser.add_argument('-l', '--log', dest='logLevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO', help='Set logging level.')
    parser.add_argument('-s', '--show', help='Show display.', action='store_true')
    parser.add_argument('-w', '--ods', help='Write output to ODS file.', action='store_true')
    parser.add_argument('-o', '--outdir', type=str, nargs='?', default=os.getcwd(), help='Output directory.')
    parser.add_argument('-c', '--correct', nargs='*', type=float, help='Polynomial correction parameters (a0, a1, a2).')
    
    args = parser.parse_args()

    # Argument Validation
    if args.brho is None and args.fref is None and args.kenergy is None and args.gamma is None:
        parser.error('You must provide one reference parameter: -f, -b, -ke, or -gam.')

    # Logging Setup
    if args.logLevel: 
        log.basicConfig(level=log.getLevelName(args.logLevel))
    
    # Handle alphap conversion
    if args.alphap and args.alphap > 1: 
        args.alphap = 1 / args.alphap**2

    print(f'Running {scriptname}...')
    log.info(f'Processing {args.datafile} for reference {args.refion}.')

    # Handle file lists vs single files
    files_to_process = []
    if 'txt' in args.datafile[0] and len(args.datafile) == 1:
        files_to_process = read_masterfile(args.datafile[0])
    else:
        files_to_process = args.datafile

    # Initialize Qt Application ONCE (Crucial for loops)
    app = None
    if args.show:
        app = QApplication.instance()
        if not app:
            app = QApplication(sys.argv)

    # Run Controller for each file
    for file in files_to_process:
        run_controller(
            data_file=file, 
            particles_to_simulate=args.filep, 
            alphap=args.alphap, 
            ref_ion=args.refion, 
            harmonics=args.harmonics, 
            brho=args.brho, 
            fref=args.fref, 
            ke=args.kenergy, 
            gam=args.gamma, 
            correct=args.correct, 
            ods=args.ods, 
            nions=args.nions,
            show=args.show,
            app=app
        )
    
    if args.show and app:
        sys.exit(app.exec_())

def run_controller(data_file, particles_to_simulate, alphap, ref_ion, harmonics, 
                   brho=None, fref=None, ke=None, gam=None, correct=None, 
                   ods=False, nions=None, show=True, app=None):
    """
    Unified controller: Calculates physics and launches the PyQt GUI.

    Parameters
    ----------
    data_file : str
        Path to experimental data.
    particles_to_simulate : str
        Path to LISE++ file.
    alphap : float
        Momentum compaction factor.
    ref_ion : str
        Reference ion string.
    harmonics : list
        List of harmonics to simulate.
    brho : float, optional
        Magnetic rigidity.
    fref : float, optional
        Reference frequency.
    ke : float, optional
        Kinetic energy.
    gam : float, optional
        Gamma factor.
    correct : list, optional
        Polynomial correction coefficients.
    ods : bool, optional
        Whether to export to ODS.
    nions : int, optional
        Number of top ions to display.
    show : bool, optional
        Whether to show the GUI.
    app : QApplication, optional
        The Qt application instance.
    """
    # 1. Calculations (Model)
    mydata = ImportData(ref_ion, alphap, filename=data_file)
    log.debug(f'Experimental data shape: {shape(mydata.experimental_data)}')
    
    mydata._set_particles_to_simulate_from_file(particles_to_simulate)
    mydata._calculate_moqs()
    
    # Calculate Reference Frequency
    mydata._calculate_srrf(fref=fref, brho=brho, ke=ke, gam=gam, correct=correct)
    log.debug(f'Reference Frequency: {mydata.ref_frequency}')
    
    # Simulate Data
    if not isinstance(harmonics, list):
        harmonics = [harmonics]
        
    mydata._simulated_data(harmonics=harmonics, brho=brho, mode='Frequency' if fref else 'Brho') 

    # 2. Filter Ions (Optional)
    if nions: 
        display_nions(nions, mydata.yield_data, mydata.nuclei_names, mydata.simulated_data_dict, ref_ion, harmonics)

    # 3. Logging & ODS Export
    sort_index = argsort(mydata.srrf)
    if ods: 
        write_arrays_to_ods(
            'Data_simulated_RionID', 
            'Data', 
            ['Name', 'freq', 'yield'], 
            (mydata.nuclei_names)[sort_index], 
            (mydata.srrf)[sort_index] * mydata.ref_frequency, 
            (mydata.yield_data)[sort_index] 
        )

    # 4. Visualization (View)
    if show and app:
        sa = CreatePyGUI(mydata.experimental_data, mydata.simulated_data_dict)
        sa.show()

def display_nions(nions, yield_data, nuclei_names, simulated_data_dict, ref_ion, harmonics):
    """Filters the top N ions by yield."""
    sorted_indices = argsort(yield_data)[::-1][:nions]
    ref_index = where(nuclei_names == ref_ion)[0]
    
    # Ensure reference ion is always included
    if len(ref_index) > 0 and ref_index[0] not in sorted_indices:
        sorted_indices = append(sorted_indices, ref_index)
        
    nuclei_names = nuclei_names[sorted_indices]
    
    for harmonic in harmonics: 
        name = f'{harmonic}'
        if name in simulated_data_dict:
            simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]

def read_masterfile(master_filename):
    """Reads list of filenames from a text file."""
    return [file.strip() for file in open(master_filename).readlines() if file.strip()]

if __name__ == '__main__':
    main()