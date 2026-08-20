import os
import re
import sys

import numpy as np
from numpy import array, polyval, sqrt, stack

from rionid.external.lisereader.reader import LISEreader
from rionid.io import handle_spectrumnpz_data
from rionid.masses import AMEData, Ring, get_ame_data, ionic_moq_u


class ImportData(object):
    """
    The core data model for RionID.

    This class handles the loading of experimental data, the physics calculations
    for ion revolution frequencies, and the simulation of expected spectra based
    on input parameters (LISE++ files, ring settings).

    Parameters
    ----------
    refion : str
        Reference ion string (e.g., '72Ge+35').
    alphap : float
        Momentum compaction factor of the ring.
    filename : str, optional
        Path to the experimental data file.
    reload_data : bool, optional
        If True, reloads raw data; otherwise loads from cache.
    circumference : float, optional
        Ring circumference in meters.
    highlight_ions : str or list, optional
        Ions to highlight in the plot.
    io_params : dict, optional
        Extra parameters for file I/O (e.g., NPZ keys).
    """

    def __init__(
        self,
        refion,
        alphap,
        filename=None,
        reload_data=None,
        circumference=None,
        highlight_ions=None,
        io_params=None,
    ):

        self.simulated_data_dict = {}
        self.particles_to_simulate = []
        self.moq = dict()
        self.protons = dict()
        self.total_mass = dict()
        self.yield_data = []

        self.highlight_ions = self._parse_highlight_ions(highlight_ions)
        self.alphap = alphap
        self.gammat = 1.0 / (self.alphap**0.5)

        self.ring = Ring("ESR", circumference)

        self.ref_ion = refion.strip()
        self._parse_ref_ion(refion)

        self.io_params = io_params or {}

        self.cache_file = self._get_cache_file_path(filename) if filename else None
        self.experimental_data = None

        if filename is not None:
            if reload_data:
                self._get_experimental_data(filename)
                self._save_experimental_data()
            else:
                try:
                    self._load_experimental_data()
                except (FileNotFoundError, IOError):
                    self._get_experimental_data(filename)

            # --- Data processing: log-safety clip + normalise ---
            if self.experimental_data is not None:
                freq, amp = self.experimental_data

                # Log-Safety (Clip negatives): ensure all values are > 0
                # for logarithmic plotting. This floor (1e-29) intentionally
                # differs from gui/plot.py's separate 1e-9 floor for the
                # same purpose.
                amp = np.maximum(amp, 1e-29)

                # Normalization: scale so the highest peak is 1.0
                max_val = np.max(amp)
                if max_val > 0:
                    amp = amp / max_val

                self.experimental_data = (freq, amp)

    def _parse_ref_ion(self, refion):
        # Regex to extract Mass(Digits), Element(Letters), Charge(Digits)
        # It handles both '98Zr+39' and '98Zr39+' inputs
        match = re.match(r"(\d+)([a-zA-Z]+).*?(\d+)", self.ref_ion)
        if match:
            self.ref_aa = int(match.group(1))
            self.ref_el = match.group(2)
            self.ref_charge = int(match.group(3))
            # Force standard format: 98Zr39+
            self.ref_ion = f"{self.ref_aa}{self.ref_el}{self.ref_charge}+"
        else:
            # Fallback parsing
            try:
                # Try splitting by '+' if it exists in the middle
                if "+" in refion and not refion.endswith("+"):
                    parts = refion.split("+")
                    self.ref_charge = int(parts[1])
                    self.ref_aa = int(re.split(r"(\d+)", parts[0])[1])
                else:
                    # Assume format like 98Zr39+
                    self.ref_charge = int(re.findall(r"\d+", refion)[-1])
                    self.ref_aa = int(re.findall(r"\d+", refion)[0])
            except Exception:
                print(f"Warning: Could not parse reference ion '{refion}'.")

    def _parse_highlight_ions(self, input_str):
        """Parses a comma-separated string of ions into a list."""
        if not input_str:
            return []
        if isinstance(input_str, list):
            return input_str
        return [x for x in re.split(r"[\s,]+", input_str.strip()) if x]

    def _get_cache_file_path(self, filename):
        """Generates the cache filename."""
        base, _ = os.path.splitext(filename)
        return f"{base}_cache.npz"

    def _get_experimental_data(self, filename):
        """Loads experimental data from various file formats."""
        base, file_extension = os.path.splitext(filename)
        ext = file_extension.lower()

        if ext == ".npz":
            self.experimental_data = handle_spectrumnpz_data(filename, **self.io_params)
        elif ext == ".root":
            raise ValueError(
                "ROOT files are not supported in this version. Please convert to NPZ/CSV."
            )

    def _save_experimental_data(self):
        """Caches loaded data to a compressed NPZ file."""
        if self.experimental_data is not None:
            frequency, amplitude_avg = self.experimental_data
            np.savez_compressed(self.cache_file, frequency=frequency, amplitude_avg=amplitude_avg)

    def _load_experimental_data(self):
        """Loads data from the cache file."""
        if os.path.exists(self.cache_file):
            data = np.load(self.cache_file, allow_pickle=True)
            frequency = data["frequency"]
            amplitude_avg = data["amplitude_avg"]
            self.experimental_data = (frequency, amplitude_avg)
        else:
            raise FileNotFoundError(
                "Cached data file not found. Please set reload_data to True to generate it."
            )

    def _set_particles_to_simulate_from_file(self, particles_to_simulate):
        """Parses the LISE++ output file and loads the (process-cached)
        AME table -- see rionid.masses.get_ame_data."""
        self.ame = get_ame_data()
        self.ame_data = self.ame.ame_table
        lise = LISEreader(particles_to_simulate)
        self.particles_to_simulate = lise.get_info_all()

    def _calculate_moqs(self):
        """Calculates mass-to-charge ratios for every candidate in
        self.particles_to_simulate, via an O(1) AME-table lookup per
        candidate.
        """
        self.moq = dict()
        self.total_mass = dict()

        for particle in self.particles_to_simulate:
            ion_name = f"{particle[1]}{particle[0]}{particle[4][-1]}+"
            ame_row = self.ame.lookup(particle[0], particle[1])
            if ame_row is None:
                continue
            qq = particle[4][-1]
            m_q = ionic_moq_u(ame_row, qq)
            self.moq[ion_name] = m_q
            self.total_mass[ion_name] = m_q * qq
            self.protons[ion_name] = ame_row[4]

    def _calculate_srrf(self, fref=None, brho=None, ke=None, gam=None, correct=None):
        """
        Calculates Simulated Relative Revolution Frequencies (SRRF).

        Applies the slip factor formula and optional polynomial correction.
        """
        self.ref_mass = AMEData.to_mev(self.moq[self.ref_ion] * self.ref_charge)
        self.ref_frequency = self.reference_frequency(fref, brho, ke, gam)
        self.srrf = array(
            [
                1 - self.alphap * (self.moq[name] - self.moq[self.ref_ion]) / self.moq[self.ref_ion]
                for name in self.moq
            ]
        )
        if correct:
            correction = polyval(array(correct), self.srrf * self.ref_frequency)
            self.srrf = self.srrf + correction / self.ref_frequency

    def _simulated_data(
        self, brho=None, harmonics=None, mode=None, sim_scalingfactor=None, nions=None
    ):
        """Generates the final simulation dictionary for plotting."""
        for harmonic in harmonics:
            ref_moq = self.moq[self.ref_ion]
            if (mode or "").lower() == "brho":
                self.brho = brho
                ref_frequency = self.ref_frequency * harmonic
            else:
                ref_frequency = self.ref_frequency
                self.brho = self.calculate_brho_relativistic(
                    ref_moq, ref_frequency, self.ring.circumference, harmonic
                )

        self.simulated_data_dict = {}
        moq_keys = list(self.moq.keys())

        yield_by_name = {}
        for p in self.particles_to_simulate:
            p_name = f"{int(p[1])}{p[0]}{int(p[4][-1])}+"
            if p_name not in yield_by_name:
                yield_by_name[p_name] = p[5]
        self.yield_data = [yield_by_name.get(key, 0) for key in moq_keys]

        self.nuclei_names = array(moq_keys)
        self.yield_data = np.array(self.yield_data, dtype=float)
        max_yield = np.max(self.yield_data)
        if max_yield > 0:
            self.yield_data /= max_yield

        if sim_scalingfactor:
            self.yield_data *= sim_scalingfactor

        for harmonic in harmonics:
            harmonic_freq = self.srrf * self.ref_frequency * harmonic
            arr_stack = stack((harmonic_freq, self.yield_data, self.nuclei_names), axis=1)
            self.simulated_data_dict[f"{harmonic}"] = arr_stack

    def calculate_brho_relativistic(self, moq, frequency, circumference, harmonic):
        """Calculates Magnetic Rigidity (Brho) from frequency."""
        actual_frequency = frequency / harmonic
        v = actual_frequency * circumference
        gamma = 1 / np.sqrt(1 - (v / AMEData.CC) ** 2)
        p = moq * AMEData.UU * gamma * (v / AMEData.CC)
        brho = (p / AMEData.CC) * 1e6
        return brho

    def reference_frequency(self, fref=None, brho=None, ke=None, gam=None):
        """Determines the reference frequency based on input mode."""
        if fref:
            return fref
        elif brho:
            return ImportData.calc_ref_rev_frequency(
                self.ref_mass, self.ring.circumference, brho=brho, ref_charge=self.ref_charge
            )
        elif ke:
            return ImportData.calc_ref_rev_frequency(
                self.ref_mass, self.ring.circumference, ke=ke, aa=self.ref_aa
            )
        elif gam:
            return ImportData.calc_ref_rev_frequency(
                self.ref_mass, self.ring.circumference, gam=gam
            )
        else:
            sys.exit("Error: No reference parameter provided.")

    @staticmethod
    def calc_ref_rev_frequency(
        ref_mass, ring_circumference, brho=None, ref_charge=None, ke=None, aa=None, gam=None
    ):
        """Static helper to calculate revolution frequency."""
        if brho is not None:
            gamma = ImportData.gamma_brho(brho, ref_charge, ref_mass)
        elif ke is not None:
            gamma = ImportData.gamma_ke(ke, aa, ref_mass)
        elif gam is not None:
            gamma = gam
        else:
            raise ValueError("Provide one of brho, ke, or gam to calculate frequency.")
        beta = ImportData.beta(gamma)
        return ImportData.velocity(beta) / ring_circumference

    @staticmethod
    def gamma_brho(brho, charge, mass):
        return sqrt(pow(brho * charge * AMEData.CC / (mass * 1e6), 2) + 1)

    @staticmethod
    def gamma_ke(ke, aa, ref_mass):
        return (ke * aa) / (ref_mass) + 1

    @staticmethod
    def beta(gamma):
        return sqrt(gamma**2 - 1) / gamma

    @staticmethod
    def velocity(beta):
        return AMEData.CC * beta

    @staticmethod
    def calc_revolution_frequency(velocity, ring_circumference):
        return velocity / ring_circumference
