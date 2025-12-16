import argparse
import os
import logging as log
import ezodf
import sys
from numpy import argsort, where, append, shape
from rionid import CreateGUI, ImportData, CreatePyGUI
from rionid.inouttools import *
from PyQt5.QtWidgets import QApplication

def main():
    
    scriptname = 'RionID' 
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required = True)

    # Main Arguments
    parser.add_argument('datafile', type = str, nargs = '+', help = 'Name of the input file with data.')
    parser.add_argument('-ap', '--alphap', type = float, help = 'Momentum compaction factor of the ring.')
    parser.add_argument('-r', '--refion', type = str, help = 'Reference ion with format NucleonsNameChargestate :=  AAXX+CC. Example: 72Ge+35, 1H+1, 238U+92...')
    parser.add_argument('-psim', '--filep', type = str, help = 'Read list of particles to simulate. LISE file or something else.')
    parser.add_argument('-hrm', '--harmonics', type = float, default = 1, nargs = '+', help = 'Harmonics to simulate.')

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
    if args.logLevel: log.basicConfig(level = log.getLevelName(args.logLevel))
    if args.outdir: outfilepath = os.path.join(args.outdir, '')

    # Easy way to handle alphap or gammat. If alphap is greater than 1, it is assumed that you are giving gammat. So here it is transformed to alphap = 1 / gammat^2
    if args.alphap > 1: args.alphap = 1 / args.alphap**2

    # Here We Go:
    print(f'Running {scriptname}... Lets see what we have in our ring ;-)')
    log.info(f'File {args.datafile} passed for processing the information of {args.refion}.')

    # If it is a txt file with files or just files introduced by the terminal
    if ('txt') in args.datafile[0]:
        datafile_list = read_masterfile(args.datafile[0])
        for datafile in datafile_list:
            controller(datafile[0], args.filep, args.harmonics, args.alphap, args.refion, args.ndivs, args.amplitude, args.show, brho = args.brho, fref = args.fref, ke = args.kenergy, out = args.outdir, harmonics = args.harmonics, gam = args.gamma, correct = args.correct, ods = args.ods, nions=args.nions)
    else:
        for file in args.datafile:
            controller2(file, args.filep, args.alphap, args.refion, args.ndivs, args.amplitude, args.show, brho = args.brho, fref = args.fref, ke = args.kenergy, out = args.outdir, harmonics = args.harmonics, gam = args.gamma, correct = args.correct, ods = args.ods, nions = args.nions)
    
def controller(data_file, particles_to_simulate, alphap, ref_ion, ndivs, amplitude, show, brho = None, fref = None, ke = None, out = None, harmonics = None, gam = None, correct = None, ods = False, nions = None):
    # Calculations
    mydata = ImportData(ref_ion, alphap, filename = data_file)
    log.debug(f'Experimental data (shape = {shape(mydata.experimental_data)}) = {mydata.experimental_data}')
    mydata._set_particles_to_simulate_from_file(particles_to_simulate)
    
    mydata._calculate_moqs()
    log.debug(f'moqs = {mydata.moq}')
    mydata._calculate_srrf(fref = fref, brho = brho, ke = ke, gam = gam, correct = correct)
    log.debug(f'Revolution (or meassured) frequency of {ref_ion} = {mydata.ref_frequency}')
    mydata._simulated_data(harmonics = harmonics) # -> simulated frecs

    log.debug(f'Simulated data (shape = {shape(mydata.simulated_data_dict[str(1.0)])}) = {mydata.simulated_data_dict}')
    log.debug(f'Simulation results (ordered by frequency) = ')
    sort_index = argsort(mydata.srrf)
    for i in sort_index:
        log.debug(f'{mydata.nuclei_names[i]} with simulated rev freq: {mydata.srrf[i] * mydata.ref_frequency} and yield: {mydata.yield_data[i]}')
    if ods: write_arrays_to_ods('Data_simulated_RionID', 'Data', ['Name', 'freq', 'yield'], (mydata.nuclei_names)[sort_index], (mydata.srrf)[sort_index] * mydata.ref_frequency, (mydata.yield_data)[sort_index] )
    log.info(f'Simulation performed. Now we are going to start the display.')

    # View
    # displaying specified amount of ions, sorted by yield
    if nions: display_nions(nions, mydata.yield_data, mydata.nuclei_names, mydata.simulated_data_dict, ref_ion, harmonics)

    mycanvas = CreateGUI(ref_ion, mydata.nuclei_names, ndivs, amplitude, show)
    mycanvas._view(mydata.experimental_data, mydata.simulated_data_dict, filename = data_file, out = out)

    log.debug(f'Plotted labels = {mycanvas.labels},{mycanvas.ref_ion}')
    log.info(f'Program has ended. I hope you have found what you were looking for. :)')

def display_nions(nions, yield_data, nuclei_names, simulated_data_dict, ref_ion, harmonics):
    sorted_indices = argsort(yield_data)[::-1][:nions]
    ref_index = where(nuclei_names == ref_ion)[0]
    if ref_index not in sorted_indices:
        sorted_indices = append(sorted_indices, ref_index)
    nuclei_names = nuclei_names[sorted_indices]
    
    for harmonic in harmonics: # for each harmonic
        name = f'{harmonic}'
        simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]

def controller2(data_file, particles_to_simulate, alphap, ref_ion, ndivs, amplitude, show, brho = None, fref = None, ke = None, out = None, harmonics = None, gam = None, correct = None, ods = False, nions = None):
    # Calculations
    mydata = ImportData(ref_ion, alphap, filename = data_file)
    mydata._set_particles_to_simulate_from_file(particles_to_simulate)
    
    mydata._calculate_moqs()
    mydata._calculate_srrf(fref = fref, brho = brho, ke = ke, gam = gam, correct = correct)
    mydata._simulated_data(harmonics = harmonics) # -> simulated frecs

    if nions: display_nions(nions, mydata.yield_data, mydata.nuclei_names, mydata.simulated_data_dict, ref_ion, harmonics)

    #pyView
    app = QApplication(sys.argv)
    sa = CreatePyGUI(mydata.experimental_data, mydata.simulated_data_dict)
    sa.show()
    sys.exit(app.exec_())

def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]

if __name__ == '__main__':
    main()