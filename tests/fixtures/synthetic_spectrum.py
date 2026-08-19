"""Synthetic (non-experimental) test fixtures for RionID's regression
suite. All data here is generated, not derived from any real measurement,
and is safe to commit/redistribute publicly.
"""
import numpy as np


def make_synthetic_spectrum(path, n_points=2000, seed=42):
    """Writes a small synthetic Schottky-like spectrum to `path` as a
    2-array .npz (arr_0=frequency [Hz], arr_1=amplitude [a.u.]), readable
    by rionid.io.handle_spectrumnpz_data's default keys.

    Not real experimental data: a flat noise floor plus three Gaussian
    lines on a 1 MHz synthetic band, used only to exercise the peak
    detection / plotting code paths deterministically and quickly in CI.
    """
    rng = np.random.default_rng(seed)
    freq = np.linspace(244.5e6, 245.5e6, n_points)
    amp = np.abs(1e-3 * (1 + 0.2 * rng.standard_normal(n_points)))
    for f0, height, width in [
        (244.71e6, 1.0, 800.0),
        (244.98e6, 0.6, 600.0),
        (245.22e6, 0.35, 700.0),
    ]:
        amp += height * np.exp(-0.5 * ((freq - f0) / width) ** 2)
    np.savez(path, arr_0=freq, arr_1=amp)
    return path


def build_ame_candidates(ame_table, n, zz_max=92):
    """Builds `n` candidate tuples in the exact shape
    `rionid.core.ImportData._calculate_moqs`/`_simulated_data` expect in
    `particles_to_simulate`: `(element_name, aa, zz, nn, [charge_states],
    yield)`. Built from REAL rows of a loaded AME table (real Z/N/A/element
    name; fully-stripped charge state) -- only "which N to pick" and the
    fixed yield=1.0 are synthetic, not the nuclide data itself.
    """
    candidates = []
    seen = set()
    for row in ame_table:
        nn, zz, aa, name = row[3], row[4], row[5], row[6]
        if zz < 1 or zz > zz_max:
            continue
        key = (name, aa)
        if key in seen:
            continue
        seen.add(key)
        candidates.append((name, aa, zz, nn, [zz], 1.0))
        if len(candidates) >= n:
            break
    return candidates
