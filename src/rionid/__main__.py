import argparse
import logging as log
import os
import sys

from numpy import append, argsort, where
from PyQt5.QtWidgets import QApplication

from rionid.gui.controller import import_controller
from rionid.gui.plot import CreatePyGUI
from rionid.io import write_arrays_to_ods
from rionid.session import candidate_rows

ESR_CIRCUMFERENCE_M = 108.36


def build_parser():
    """Build the parser used by the ``python -m rionid`` interface."""
    parser = argparse.ArgumentParser(description="RionID: Ring-stored ion Identification")
    modes = parser.add_mutually_exclusive_group(required=True)

    # Main Arguments
    parser.add_argument("datafile", type=str, nargs="+", help="Name of the input file with data.")
    parser.add_argument(
        "-ap",
        "--alphap",
        type=float,
        required=True,
        help="Momentum compaction factor of the ring.",
    )
    parser.add_argument(
        "-r",
        "--refion",
        type=str,
        required=True,
        help="Reference ion (format: AAXXCC+, e.g., 72Ge35+).",
    )
    parser.add_argument(
        "-psim",
        "--filep",
        type=str,
        default="",
        help="LISE++ candidate list (optional with --add-ions or --search-ion).",
    )
    parser.add_argument(
        "-hrm", "--harmonics", type=float, default=[1.0], nargs="+", help="Harmonics to simulate."
    )

    parser.add_argument("--add-ions", default="", help="Additional ions, e.g. '205Au76+,205Au75+'.")
    parser.add_argument("--search-ion", default="", help="Search AME2020 m/q near this ion.")
    parser.add_argument(
        "--moq-tolerance-percent",
        type=float,
        default=1.0,
        help="Relative m/q half-width, in percent (default: 1).",
    )
    parser.add_argument("--charge-min", type=int, default=1)
    parser.add_argument(
        "--charge-max", type=int, default=None, help="Maximum q; default is each nuclide's Z."
    )
    parser.add_argument(
        "--display-filter", default="", help="Show only an ion, isotope or element."
    )
    parser.add_argument("--highlight-ions", default="", help="Ion, isotope or element selectors.")
    parser.add_argument(
        "--show-uncertainties", action="store_true", help="Show ±1σ mass-input bands."
    )
    parser.add_argument(
        "--frequency-scale-hz-per-100kev",
        type=float,
        help="Optional empirical conversion X Hz per 100 keV.",
    )
    parser.add_argument(
        "--uncertainty-harmonic",
        type=float,
        default=1.0,
        help="Harmonic/multiplier at which the empirical scale applies.",
    )
    parser.add_argument("--frequency-key", default="arr_0", help="Frequency array key in the NPZ.")
    parser.add_argument("--amplitude-key", default="arr_1", help="Amplitude array key in the NPZ.")

    # Secondary Arguments
    parser.add_argument(
        "-n", "--nions", type=int, help="Number of ions to display, sorted by yield."
    )

    # Arguments for Each Mode (Exclusive)
    modes.add_argument("-b", "--brho", type=float, help="Brho value [Tm] (Isochronous mode).")
    modes.add_argument(
        "-ke", "--kenergy", type=float, help="Kinetic energy [MeV/u] (Isochronous mode)."
    )
    modes.add_argument("-gam", "--gamma", type=float, help="Lorentz factor gamma.")
    modes.add_argument(
        "-f", "--fref", type=float, help="Revolution frequency [Hz] (Standard mode)."
    )

    # Visualization & Output
    parser.add_argument(
        "-d", "--ndivs", type=int, default=4, help="Number of divisions (Deprecated)."
    )
    parser.add_argument(
        "-am", "--amplitude", type=int, default=0, help="Display options (0=const, else=scaled)."
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="logLevel",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set logging level.",
    )
    parser.add_argument("-s", "--show", help="Show display.", action="store_true")
    parser.add_argument("-w", "--ods", help="Write output to ODS file.", action="store_true")
    parser.add_argument(
        "-o", "--outdir", type=str, nargs="?", default=os.getcwd(), help="Output directory."
    )
    parser.add_argument(
        "-c",
        "--correct",
        nargs="*",
        type=float,
        help="Polynomial correction parameters (a0, a1, a2).",
    )

    parser.add_argument(
        "--circumference",
        type=float,
        default=ESR_CIRCUMFERENCE_M,
        metavar="METRES",
        help=f"Ring circumference in metres (default: ESR, {ESR_CIRCUMFERENCE_M}).",
    )
    return parser


def main():
    """
    Main entry point for the command-line interface (CLI).

    Parses arguments, initializes logging, and dispatches the analysis controller
    for one or more data files.
    """
    scriptname = "RionID"
    parser = build_parser()

    args = parser.parse_args()

    if not (args.filep or args.add_ions or args.search_ion):
        parser.error("Provide --filep, --add-ions or --search-ion.")
    # Argument Validation
    if args.circumference <= 0:
        parser.error("--circumference must be greater than zero.")

    # Logging Setup
    if args.logLevel:
        log.basicConfig(level=log.getLevelName(args.logLevel))

    # Handle alphap conversion
    if args.alphap and args.alphap > 1:
        args.alphap = 1 / args.alphap**2

    print(f"Running {scriptname}...")
    log.info(f"Processing {args.datafile} for reference {args.refion}.")

    # Handle file lists vs single files
    files_to_process = []
    if "txt" in args.datafile[0] and len(args.datafile) == 1:
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
            circumference=args.circumference,
            brho=args.brho,
            fref=args.fref,
            ke=args.kenergy,
            gam=args.gamma,
            correct=args.correct,
            ods=args.ods,
            nions=args.nions,
            show=args.show,
            app=app,
            manual_ions=args.add_ions,
            search_ion=args.search_ion,
            moq_tolerance_percent=args.moq_tolerance_percent,
            charge_min=args.charge_min,
            charge_max=args.charge_max,
            display_filter=args.display_filter,
            highlight_ions=args.highlight_ions,
            show_uncertainties=args.show_uncertainties,
            frequency_scale_hz_per_100kev=args.frequency_scale_hz_per_100kev,
            uncertainty_harmonic=args.uncertainty_harmonic,
            io_params={"frequency_key": args.frequency_key, "amplitude_key": args.amplitude_key},
            output_dir=args.outdir,
        )

    if args.show and app:
        sys.exit(app.exec_())


def run_controller(
    data_file,
    particles_to_simulate,
    alphap,
    ref_ion,
    harmonics,
    circumference=ESR_CIRCUMFERENCE_M,
    brho=None,
    fref=None,
    ke=None,
    gam=None,
    correct=None,
    ods=False,
    nions=None,
    show=True,
    app=None,
    manual_ions="",
    search_ion="",
    moq_tolerance_percent=1.0,
    charge_min=1,
    charge_max=None,
    display_filter="",
    highlight_ions="",
    show_uncertainties=False,
    frequency_scale_hz_per_100kev=None,
    uncertainty_harmonic=1.0,
    io_params=None,
    output_dir=None,
):
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
    circumference : float, optional
        Ring circumference in metres. Defaults to the ESR value of 108.36 m.
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
    if fref is not None:
        mode, value = "Frequency", fref
    elif brho is not None:
        mode, value = "Brho", brho
    elif ke is not None:
        mode, value = "Kinetic Energy", ke
    else:
        mode, value = "Gamma", gam
    if not isinstance(harmonics, (list, tuple)):
        harmonics = [harmonics]
    mydata = import_controller(
        datafile=data_file,
        filep=particles_to_simulate or "",
        alphap=alphap,
        refion=ref_ion,
        harmonics=list(harmonics),
        circumference=circumference,
        mode=mode,
        value=value,
        correct=correct,
        nions=nions,
        reload_data=True,
        io_params=io_params or {"frequency_key": "arr_0", "amplitude_key": "arr_1"},
        output_path=None,
        write_cache=False,
        manual_ions=manual_ions,
        search_ion=search_ion,
        moq_tolerance_percent=moq_tolerance_percent,
        charge_min=charge_min,
        charge_max=charge_max,
        display_filter=display_filter,
        highlight_ions=highlight_ions,
        show_uncertainties=show_uncertainties,
        frequency_scale_hz_per_100kev=frequency_scale_hz_per_100kev,
        uncertainty_harmonic=uncertainty_harmonic,
    )
    if ods:
        from pathlib import Path

        destination = Path(output_dir or os.getcwd())
        destination.mkdir(parents=True, exist_ok=True)
        rows = candidate_rows(mydata)
        columns = list(rows[0])
        write_arrays_to_ods(
            str(destination / "Data_simulated_RionID.ods"),
            "Candidates",
            columns,
            *[["" if row[column] is None else row[column] for row in rows] for column in columns],
        )
    if show and app:
        window = CreatePyGUI()
        window.updateData(mydata)
        window.show()
        if not hasattr(app, "rionid_windows"):
            app.rionid_windows = []
        app.rionid_windows.append(window)

    return mydata


def display_nions(nions, yield_data, nuclei_names, simulated_data_dict, ref_ion, harmonics):
    """Filters the top N ions by yield."""
    sorted_indices = argsort(yield_data)[::-1][:nions]
    ref_index = where(nuclei_names == ref_ion)[0]

    # Ensure reference ion is always included
    if len(ref_index) > 0 and ref_index[0] not in sorted_indices:
        sorted_indices = append(sorted_indices, ref_index)

    nuclei_names = nuclei_names[sorted_indices]

    for harmonic in harmonics:
        name = f"{harmonic}"
        if name in simulated_data_dict:
            simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]


def read_masterfile(master_filename):
    """Reads list of filenames from a text file."""
    return [file.strip() for file in open(master_filename).readlines() if file.strip()]


if __name__ == "__main__":
    main()
