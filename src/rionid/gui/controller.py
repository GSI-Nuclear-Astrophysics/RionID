import numpy as np
from loguru import logger
from numpy import append, argsort, where

from rionid.candidates import parse_ion, split_ions
from rionid.core import ImportData
from rionid.masses import AMEData


def import_controller(
    datafile=None,
    filep=None,
    alphap=None,
    refion=None,
    harmonics=None,
    nions=None,
    amplitude=None,
    circumference=None,
    mode=None,
    value=None,
    reload_data=None,
    highlight_ions=None,
    io_params=None,
    sim_scalingfactor=None,
    correct=None,
    output_path="simulation_result.out",
    write_cache=True,
    ionic_inputs=None,
    manual_ions="",
    search_ion="",
    moq_tolerance_percent=1.0,
    charge_min=1,
    charge_max=None,
    display_filter="",
    frequency_scale_hz_per_100kev=None,
    uncertainty_harmonic=1.0,
    show_uncertainties=False,
):
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
    highlight_ions : str, optional
        Comma-separated string of ion names to highlight in the plot (e.g., '72Ge+35, 74Se+34').
    io_params : dict, optional
        Dictionary containing specific I/O parameters (e.g., keys for NPZ files,
        histogram names for ROOT files).
    sim_scalingfactor : float, optional
        Factor to scale the simulated yield/amplitude.
    correct : list of float, optional
        Coefficients [a0, a1, a2] for a second-order polynomial correction
        applied to the simulated frequencies.
    output_path : str or None, optional
        Text output filename; None disables automatic export.
    write_cache : bool, optional
        Whether reloading writes a spectrum cache beside the source.
    ionic_inputs : dict or None, optional
        Resolved ionic masses and charges from a completed run snapshot;
        None uses the current AME2020 table and binding-energy correction.

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
        mode = "Brho" if mode == "Bρ" else mode
        if float(alphap) > 1:
            alphap = 1 / float(alphap) ** 2  # handling alphap and gammat
        fref = brho = ke = gam = None
        if mode == "Frequency":
            fref = float(value)
        elif mode == "Brho":
            brho = float(value)
        elif mode == "Kinetic Energy":
            ke = float(value)
        elif mode == "Gamma":
            gam = float(value)
        # Calculations
        mydata = ImportData(
            refion,
            float(alphap),
            filename=datafile,
            reload_data=reload_data,
            circumference=circumference,
            highlight_ions=highlight_ions,
            io_params=io_params,
            write_cache=write_cache,
        )

        mydata.display_filter = display_filter
        mydata.show_uncertainties = show_uncertainties
        mydata.frequency_scale_hz_per_100kev = frequency_scale_hz_per_100kev
        mydata.uncertainty_harmonic = uncertainty_harmonic
        if ionic_inputs is None:
            if filep:
                mydata._set_particles_to_simulate_from_file(filep)
                mydata._calculate_moqs()
            if manual_ions or search_ion or not filep:
                mydata._extend_candidates(
                    manual_ions, search_ion, moq_tolerance_percent, charge_min, charge_max
                )
            if mydata.ref_ion not in mydata.moq:
                mydata._extend_candidates()
        else:
            if filep:
                mydata._set_particles_to_simulate_from_file(filep, load_masses=False)
            allowed_names = {f"{p[1]}{p[0]}{p[4][-1]}+" for p in mydata.particles_to_simulate}
            allowed_names.update(parse_ion(name).name for name in split_ions(manual_ions))
            allowed_names.add(mydata.ref_ion)
            if not search_ion and not set(ionic_inputs) <= allowed_names:
                raise ValueError("Saved ionic inputs do not match the candidate list")
            for ion, entry in ionic_inputs.items():
                parsed = parse_ion(ion)
                moq, mass = float(entry["moq_u_per_e"]), float(entry["mass_u"])
                if not np.isfinite(moq) or not np.isfinite(mass) or moq <= 0 or mass <= 0:
                    raise ValueError(f"Invalid saved ionic mass for {ion}")
                if int(entry["Z"]) != parsed.Z or not np.isclose(mass, moq * parsed.q, rtol=1e-12):
                    raise ValueError(f"Inconsistent saved ionic coordinate for {ion}")
                metadata = entry.get("metadata", {})
                sigma = metadata.get("mass_uncertainty_keV")
                if sigma is not None and (not np.isfinite(sigma) or sigma < 0):
                    raise ValueError(f"Invalid saved mass uncertainty for {ion}")
                mydata.moq[ion] = moq
                mydata.total_mass[ion] = mass
                mydata.protons[ion] = int(entry["Z"])
                mydata.candidate_metadata[ion] = dict(
                    entry.get(
                        "metadata",
                        {
                            "source": "lise",
                            "sources": ["lise"],
                            "mass_uncertainty_keV": None,
                            "mass_extrapolated": None,
                            "relative_yield_known": True,
                        },
                    )
                )
        mydata._calculate_srrf(fref=fref, brho=brho, ke=ke, gam=gam, correct=correct)

        if isinstance(harmonics, str):
            harmonics = [float(h.strip()) for h in harmonics.replace(",", " ").split()]
        elif not isinstance(harmonics, list):
            harmonics = [float(harmonics)]

        mydata._simulated_data(
            brho=brho,
            harmonics=harmonics,
            mode=mode,
            sim_scalingfactor=sim_scalingfactor,
            nions=nions,
        )
        # "Outputs"
        if nions:
            display_nions(
                int(nions),
                mydata.yield_data,
                mydata.nuclei_names,
                mydata.simulated_data_dict,
                mydata.ref_ion,
                harmonics,
            )

        if nions:
            for harmonic, full_rows in mydata.all_simulated_data_dict.items():
                visible = set(mydata.simulated_data_dict[harmonic][:, 2])
                for row in full_rows:
                    sources = mydata.candidate_metadata.get(str(row[2]), {}).get("sources", [])
                    if any(source in sources for source in ("manual", "ame_search", "reference")):
                        visible.add(str(row[2]))
                mydata.simulated_data_dict[harmonic] = full_rows[
                    [str(row[2]) in visible for row in full_rows]
                ]

        if output_path is not None:
            logger.info(f"Saving simulation results to {output_path}")
            sort_index = argsort(mydata.srrf)
            save_simulation_results(mydata, mode, harmonics, sort_index, filename=output_path)

        return mydata  # Returns the simulated spectrum data

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
        name = f"{harmonic}"
        if name in simulated_data_dict:
            simulated_data_dict[name] = simulated_data_dict[name][sorted_indices]


def save_simulation_results(mydata, mode, harmonics, sort_index, filename="simulation_result.out"):
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
    with open(filename, "w") as file:
        # Writing harmonics and brho information
        brho = mydata.brho
        for harmonic in harmonics:
            header0 = f"Harmonic: {harmonic} , Bp: {brho:.6f} [Tm]"
            logger.info(header0)
            file.write(header0 + "\n")

        # Writing the header for the data table
        header1 = f"{'ion':<15}{'fre[Hz]':<30}{'yield [pps]':<15}{'m/q [u]':<15}{'m [eV]':<15}"
        file.write(header1 + "\n")
        file.write("-" * len(header1) + "\n")
        logger.info(header1)

        # Writing the sorted simulation results
        for i in sort_index:
            ion = mydata.nuclei_names[i]
            fre = mydata.srrf[i] * mydata.ref_frequency
            if mode in ("Brho", "Bρ"):
                fre *= harmonics[-1]
            yield_ = mydata.yield_data[i]
            moq = mydata.moq[ion]
            mass_u = mydata.total_mass[ion]
            mass = AMEData.to_mev(mass_u) * 1e6
            result_line = f"{ion:<15}{fre:<30.10f}{yield_:<15.4e}{moq:<15.12f}{mass:<15.3f}"
            logger.info(result_line)
            file.write(result_line + "\n")
