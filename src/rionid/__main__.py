import argparse
import os
import logging as log
import sys
from numpy import argsort, where, append, shape
from PyQt5.QtWidgets import QApplication

from rionid import ImportData, CreatePyGUI
from rionid.io import write_arrays_to_ods

def main():
    
    scriptname = 'RionID' 
    parser = argparse.ArgumentParser(description="RionID: Ring-stored ion Identification")
    modes = parser.add_mutually_exclusive_group(required = True)

    # Main Arguments
    parser.add_argument('datafile', type = str, nargs = '+', help = 'Name of the input file with data.')
    parser.add_argument('-ap', '--alphap', type = float, help = 'Momentum compaction factor of the ring.')
    parser.add_argument('-r', '--refion', type = str, help = 'Reference ion with format NucleonsNameChargestate :=  AAXX+CC. Example: 72Ge+35, 1H+1, 238U+92...')
    parser.add_argument('-psim', '--filep', type = str, help = 'Read list of particles to simulate. LISE file or something else.')
    parser.add_argument('-hrm', '--harmonics', type = float, default = [1], nargs = '+', help = 'Harmonics to simulate.')

    # Secondary Arguments
    parser.add_argument('-n', '--nions', type = int, help = 'Number of ions to display, sorted by yield (highest)')

    # Arguments for Each Mode (Exclusive)
    modes.add_argument('-b', '--brho', type = float, help = 'Brho value of the reference nucleus at ESR (isochronous mode).')
    modes.add_argument('-ke', '--kenergy', type = float, help = 'Kinetic energy of reference nucleus at ESR (isochronous mode).')
    modes.add_argument('-gam', '--gamma', type = float, help = 'Lorentz factor gamma of the reference particle')
    modes.add_argument('-f', '--fref', type = float, help = 'Revolution frequency of the reference particle (standard mode).')
    
    # Arguments for the Visualization
    parser.add_argument('-d', '--ndivs', type = int, default = 4, help = 'Number of divisions in the display.')
    parser.add_argument('-am', '--amplitude', type = int, default = 0, help = 'Display of srf data options. 0 -> constant height, else->scaled.')
    
    # Actions
    parser.add_argument('-l', '--log', dest = 'logLevel', choices = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default = 'INFO', help = 'Set the logging level.')
    parser.add_argument('-s', '--show', help = 'Show display. If not, save root file and close display', action = 'store_true')
    parser.add_argument('-w', '--ods', help = 'Write ods.', action = 'store_true')

    parser.add_argument('-o', '--outdir', type = str, nargs = '?', default = os.getcwd(), help = 'Output directory.')
    parser.add_argument('-c', '--correct', nargs = '*', type = float, help = 'Correct simulated spectrum following a polynomial fit with paremeters given here')
    
    args = parser.parse_args()

    # Checking for Argument Errors
    if args.brho is None and args.fref is None and args.kenergy is None and args.gamma is None:
        parser.error('Please introduce the revolution frequency of the reference nucleus or the brho parameter or ke/aa or gamma.')

    # Extra Details
    if args.logLevel: 
        log.basicConfig(level = log.getLevelName(args.logLevel))
    if args.outdir: 
        outfilepath = os.path.join(args.outdir, '')

    # Easy way to handle alphap or gammat. If alphap is greater than 1, it is assumed that you are giving gammat. So here it is transformed to alphap = 1 / gammat^2
    if args.alphap > 1: 
        args.alphap = 1 / args.alphap**2

    # Here We Go:
    print(f'Running {scriptname}... Lets see what we have in our ring ;-)')
    log.info(f'File {args.datafile} passed for processing the information of {args.refion}.')

    # If it is a txt file with files or just files introduced by the terminal
    files_to_process = []
    if 'txt' in args.datafile[0] and len(args.datafile) == 1:
        files_to_process = read_masterfile(args.datafile[0])
    else:
        files_to_process = args.datafile

    # Initialize Qt Application ONCE
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
    

def run_controller(data_file, particles_to_simulate, alphap, ref_ion, ndivs, amplitude, show, brho = None, fref = None, ke = None, out = None, harmonics = None, gam = None, correct = None, ods = False, nions = None):
    # Calculations
    mydata = ImportData(ref_ion, alphap, filename = data_file)
    log.debug(f'Experimental data shape: {shape(mydata.experimental_data)}')

    mydata._set_particles_to_simulate_from_file(particles_to_simulate)
    mydata._calculate_moqs()

    mydata._calculate_srrf(fref = fref, brho = brho, ke = ke, gam = gam, correct = correct)
    log.debug(f'Reference Frequency: {mydata.ref_frequency}')
    if not isinstance(harmonics, list):
        harmonics = [harmonics]
        
    mydata._simulated_data(harmonics=harmonics, brho=brho, mode='Frequency' if fref else 'Brho')

    if nions: 
        display_nions(nions, mydata.yield_data, mydata.nuclei_names, mydata.simulated_data_dict, ref_ion, harmonics)

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
        app.exec_()

def display_nions(nions, yield_data, nuclei_names, simulated_data_dict, ref_ion, harmonics):
    """Filters the top N ions by yield."""
    sorted_indices = argsort(yield_data)[::-1][:nions]
    ref_index = where(nuclei_names == ref_ion)[0]
    if ref_index not in sorted_indices:
        sorted_indices = append(sorted_indices, ref_index)
    nuclei_names = nuclei_names[sorted_indices]
    
    for harmonic in harmonics: # for each harmonic
        name = f'{harmonic}'
        if name in simulated_data_dict:
            simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]

def read_masterfile(master_filename):
    # Reads list of filenames from a text file.
    return [file.strip() for file in open(master_filename).readlines() if file.strip()]

if __name__ == '__main__':
    main()