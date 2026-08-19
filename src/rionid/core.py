from numpy import polyval, array, stack, append, sqrt
import sys
import re
import os
import numpy as np
from scipy.signal import find_peaks, peak_widths

from rionid.masses import Ring, AMEData, get_ame_data, ionic_moq_u
from rionid.external.lisereader.reader import LISEreader

from rionid.io import (
    read_psdata,
    handle_read_tdsm_bin,
    handle_spectrumnpz_data,
)


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
    peak_threshold_pct : float, optional
        Peak detection threshold (0.0-1.0).
    min_distance : float, optional
        Minimum distance between peaks.
    matching_freq_min : float, optional
        Minimum frequency for peak matching.
    matching_freq_max : float, optional
        Maximum frequency for peak matching.
    io_params : dict, optional
        Extra parameters for file I/O (e.g., NPZ keys).
    """

    def __init__(self, refion, alphap, filename=None, reload_data=None, circumference=None,
                 highlight_ions=None, peak_threshold_pct=0.05, min_distance=10,
                 matching_freq_min=None, matching_freq_max=None, io_params=None):

        self.simulated_data_dict = {}
        self.particles_to_simulate = []
        self.moq = dict()
        self.protons = dict()
        self.total_mass = dict()
        self.yield_data = []
        
        self.highlight_ions = self._parse_highlight_ions(highlight_ions)
        self.alphap = alphap
        self.gammat = 1.0 / (self.alphap ** 0.5)

        self.ring = Ring('ESR', circumference)

        self.ref_ion = refion.strip()
        self._parse_ref_ion(refion)
        
        # Physics / Matching Params
        self.peak_threshold_pct = float(peak_threshold_pct) if peak_threshold_pct else 0.05
        self.min_distance = float(min_distance) if min_distance else 10
        self.matching_freq_min = matching_freq_min
        self.matching_freq_max = matching_freq_max
        self.io_params = io_params or {}

        # Results containers
        self.peak_freqs = []
        self.peak_heights = []

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
                # for logarithmic plotting. NOTE: the floor here (1e-29)
                # intentionally does not match gui/plot.py's separate
                # 1e-9 floor for the same purpose -- pre-existing,
                # unchanged, logged in docs/LEGACY_BEHAVIOUR.md, not
                # fixed as part of this removal.
                amp = np.maximum(amp, 1e-29)

                # Normalization: scale so the highest peak is 1.0
                max_val = np.max(amp)
                if max_val > 0:
                    amp = amp / max_val

                self.experimental_data = (freq, amp)

            # Process peaks after loading and processing
            self.detect_peaks_and_widths()
    
    def _parse_ref_ion(self, refion):
        # Regex to extract Mass(Digits), Element(Letters), Charge(Digits)
        # It handles both '98Zr+39' and '98Zr39+' inputs
        match = re.match(r'(\d+)([a-zA-Z]+).*?(\d+)', self.ref_ion)
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
                if '+' in refion and not refion.endswith('+'):
                    parts = refion.split('+')
                    self.ref_charge = int(parts[1])
                    self.ref_aa = int(re.split(r'(\d+)', parts[0])[1])
                else:
                    # Assume format like 98Zr39+
                    self.ref_charge = int(re.findall(r'\d+', refion)[-1])
                    self.ref_aa = int(re.findall(r'\d+', refion)[0])
            except:
                print(f"Warning: Could not parse reference ion '{refion}'.")

    def _parse_highlight_ions(self, input_str):
        """Parses a comma-separated string of ions into a list."""
        if not input_str: return []
        if isinstance(input_str, list): return input_str
        return [x.strip() for x in input_str.split(',') if x.strip()]

    def _get_cache_file_path(self, filename):
        """Generates the cache filename."""
        base, _ = os.path.splitext(filename)
        return f"{base}_cache.npz"
    
    def _get_experimental_data(self, filename):
        """Loads experimental data from various file formats."""
        base, file_extension = os.path.splitext(filename)
        ext = file_extension.lower()

        if ext == '.csv':
            self.experimental_data = read_psdata(filename, dbm=False)
        elif ext in ['.bin_fre', '.bin_time', '.bin_amp']:
            self.experimental_data = handle_read_tdsm_bin(filename)
        elif ext == '.npz':
            self.experimental_data = handle_spectrumnpz_data(filename, **self.io_params)
        elif ext == '.root':
            raise ValueError("ROOT files are not supported in this version. Please convert to NPZ/CSV.")

    def detect_peaks_and_widths(self):
        """
        Detects peaks in the experimental spectrum using scipy.signal.find_peaks.
        
        Updates `self.peak_freqs` and `self.peak_heights`.
        """
        if self.experimental_data is None: return
        freq, amp = self.experimental_data
        
        rel_height = max(0.0, min(self.peak_threshold_pct, 1.0))
        height_thresh = np.max(amp) * rel_height
        
        peaks, _ = find_peaks(
            amp,
            height=height_thresh,
            distance=self.min_distance,
            prominence=height_thresh * 0.2,
            width=1
        )
        
        # Filter by frequency window
        peak_freqs = freq[peaks]
        mask = np.ones_like(peaks, dtype=bool)
        if self.matching_freq_min is not None:
            mask &= (peak_freqs >= self.matching_freq_min)
        if self.matching_freq_max is not None:
            mask &= (peak_freqs <= self.matching_freq_max)
            
        self.peak_freqs = peak_freqs[mask]
        self.peak_heights = amp[peaks][mask]
        self.peak_widths_freq = np.zeros_like(self.peak_freqs) 

    def _save_experimental_data(self):
        """Caches loaded data to a compressed NPZ file."""
        if self.experimental_data is not None:
            frequency, amplitude_avg = self.experimental_data
            np.savez_compressed(self.cache_file, frequency=frequency, amplitude_avg=amplitude_avg)                        
            
    def _load_experimental_data(self):
        """Loads data from the cache file."""
        if os.path.exists(self.cache_file):
            data = np.load(self.cache_file, allow_pickle=True)
            frequency = data['frequency']
            amplitude_avg = data['amplitude_avg']
            self.experimental_data = (frequency, amplitude_avg)
        else:
            raise FileNotFoundError("Cached data file not found. Please set reload_data to True to generate it.")

    def _set_particles_to_simulate_from_file(self, particles_to_simulate):
        """Parses the LISE++ output file and loads the (process-cached)
        AME table -- see rionid.masses.get_ame_data."""
        self.ame = get_ame_data()
        self.ame_data = self.ame.ame_table
        lise = LISEreader(particles_to_simulate)
        self.particles_to_simulate = lise.get_info_all()

    def _calculate_moqs(self):
        """Calculates mass-to-charge ratios for every candidate in
        self.particles_to_simulate, via an O(1) AME-table lookup (see
        docs/PERFORMANCE_BASELINE.md; previously an
        O(N x AME-table-size) linear scan per candidate). The `particles=`
        parameter this method used to accept was dead code -- confirmed
        by repo-wide grep that every call site (__main__.py,
        gui/controller.py, gui/inputs.py) always calls it with no
        arguments -- and has been removed.
        """
        self.moq = dict()
        self.total_mass = dict()

        for particle in self.particles_to_simulate:
            ion_name = f'{particle[1]}{particle[0]}{particle[4][-1]}+'
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
        self.srrf = array([1 - self.alphap * (self.moq[name] - self.moq[self.ref_ion]) / self.moq[self.ref_ion]
                           for name in self.moq])
        if correct:
            correction = polyval(array(correct), self.srrf * self.ref_frequency)
            self.srrf = self.srrf + correction / self.ref_frequency

    def _simulated_data(self, brho=None, harmonics=None, mode=None, sim_scalingfactor=None, nions=None):
        """Generates the final simulation dictionary for plotting."""
        for harmonic in harmonics:
            ref_moq = self.moq[self.ref_ion]
            if mode == 'brho':
                self.brho = brho
                ref_frequency = self.ref_frequency * harmonic
            else:
                ref_frequency = self.ref_frequency
                self.brho = self.calculate_brho_relativistic(ref_moq, ref_frequency, self.ring.circumference, harmonic)

        self.simulated_data_dict = {}
        self.yield_data = []
        moq_keys = list(self.moq.keys())
        
        for key in moq_keys:
            found = False
            for p in self.particles_to_simulate:
                p_name = f"{int(p[1])}{p[0]}{int(p[4][-1])}+"
                if p_name == key:
                    self.yield_data.append(p[5])
                    found = True
                    break
            if not found: self.yield_data.append(0)

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
            self.simulated_data_dict[f'{harmonic}'] = arr_stack
    
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
        if fref: return fref
        elif brho: return ImportData.calc_ref_rev_frequency(self.ref_mass, self.ring.circumference, brho=brho, ref_charge=self.ref_charge)
        elif ke: return ImportData.calc_ref_rev_frequency(self.ref_mass, self.ring.circumference, ke=ke, aa=self.ref_aa)
        elif gam: return ImportData.calc_ref_rev_frequency(self.ref_mass, self.ring.circumference, gam=gam)
        else: sys.exit('Error: No reference parameter provided.')
            
    @staticmethod
    def calc_ref_rev_frequency(ref_mass, ring_circumference, brho=None, ref_charge=None, ke=None, aa=None, gam=None):
        """Static helper to calculate revolution frequency."""
        if brho: gamma = ImportData.gamma_brho(brho, ref_charge, ref_mass)
        elif ke: gamma = ImportData.gamma_ke(ke, aa, ref_mass)
        elif gam: gamma = gam
        beta = ImportData.beta(gamma)
        return ImportData.velocity(beta) / ring_circumference

    @staticmethod
    def gamma_brho(brho, charge, mass): return sqrt(pow(brho * charge * AMEData.CC / (mass * 1e6), 2)+1)
    @staticmethod
    def gamma_ke(ke, aa, ref_mass): return (ke * aa) / (ref_mass) + 1
    @staticmethod
    def beta(gamma): return sqrt(gamma**2 - 1) / gamma
    @staticmethod
    def velocity(beta): return AMEData.CC * beta
    @staticmethod
    def calc_revolution_frequency(velocity, ring_circumference): return velocity / ring_circumference