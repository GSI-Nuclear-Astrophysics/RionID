# RionID (Ring-stored ion IDentification)

[![Documentation](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://GSI-Nuclear-Astrophysics.github.io/RionID)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8169341.svg)](https://doi.org/10.5281/zenodo.8169341)
[![PyPI version](https://badge.fury.io/py/rionid.svg)](https://badge.fury.io/py/rionid)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**RionID** is a Python GUI application for identifying stored highly
charged ions in storage-ring Schottky spectra. Given a reference ion, a
momentum-compaction factor, and a candidate ion list, it computes expected
revolution-frequency patterns, applies an empirical polynomial correction
anchored to known lines, and overlays the result on the experimental
spectrum for expert assessment.

<div align="center">
  <img src="https://github.com/GSI-Nuclear-Astrophysics/rionid/raw/master/docs/img/rionid.png?raw=true" width="50%">
</div>

## Scope and non-goals

RionID performs **deterministic forward simulation and expert-guided
overlay** — it computes and displays where candidate ions are expected to
appear, and lets you compare that against the data. **It does not do
automatic or autonomous ion identification.** There is no automatic
peak-to-ion assignment, no classifier, and no hidden ranking logic — every
displayed candidate is one you asked for, and every accepted assignment is
one you make yourself.

## Features

*   **Pure Python:** no ROOT dependencies.
*   **Reference-anchored forward simulation:** candidate revolution
    frequencies from ionic mass-to-charge ratios and a user-specified
    reference ion, projected to any requested harmonics.
*   **Polynomial residual correction:** an explicit, user-supplied
    quadratic correction in revolution-frequency space, applied
    consistently across harmonics.
*   **Interactive spectrum overlay:** pan, zoom, and inspect candidate ion
    labels against 1D experimental spectra.
*   **Standalone:** bundles `lisereader` (GPL-3.0) for LISE++ candidate-list
    import without extra dependency management.

## Installation

### From PyPI (recommended)

```bash
pip install rionid
```
PyPI serves the latest release.

### From source (development)

```bash
git clone https://github.com/GSI-Nuclear-Astrophysics/rionid.git
cd rionid
pip install -e ".[dev]"
```

The GUI stores local settings and file paths in `parameters_cache.toml`.
That file is ignored by Git. A sanitized template is provided as
`parameters_cache.example.toml`; copy it locally if you want to reset or
prepare the settings manually.

## Quick start

```bash
rionid
```
Fill in a reference ion, momentum-compaction factor, exactly one
reference-frequency value (frequency, Brho, kinetic energy, or gamma),
a candidate list (LISE++ `.lpp` output), and a spectrum file, then run the
simulation from the window. The GUI uses the ESR circumference of 108.36 m.

`datafile.npz` needs `arr_0`/`arr_1` keys (frequency, amplitude) by
default, or any two array keys mapped via the GUI's key-selection dialog.

For command-line use, run `rionid-cli --help` or `python3 -m rionid
--help`. The LISE++ candidate list, reference ion, momentum-compaction
factor, and exactly one reference-frequency mode are required. The ring
circumference defaults to the ESR value of 108.36 m and can be changed
with `--circumference`.

A small redistributable candidate list is included for testing the full
CLI workflow:

```bash
python examples/quickstart.py
rionid-cli examples/data/quickstart_synthetic.npz \
  -psim examples/data/candidates.lpp \
  -r 72Ge32+ -ap 0.189 -f 1930000 -hrm 127
```

## Parameter reference

| Parameter | CLI flag | Meaning |
|---|---|---|
| Data file | positional | Spectrum file (`.npz`, `.csv`, `.bin_fre`/`.bin_time`/`.bin_amp`) |
| Reference ion | `-r`, `--refion` | e.g. `72Ge+32` — sets the frequency-model anchor |
| Momentum compaction | `-ap`, `--alphap` | α_p; values `>1` are treated as γ_t and converted (`α_p = 1/γ_t²`) |
| Candidate list | `-psim`, `--filep` | LISE++ output file |
| Ring circumference | `--circumference` | Metres; defaults to ESR (108.36 m) |
| Harmonics | `-hrm`, `--harmonics` | One or more harmonic orders to display |
| Reference frequency mode | `-f`/`-b`/`-ke`/`-gam` | Exactly one of: frequency [Hz], Brho [Tm], kinetic energy [MeV/u], Lorentz γ |
| Polynomial correction | `-c`, `--correct` | `A B C` coefficients (quadratic, linear, constant), Hz-based, `numpy.polyval` order |
| Top-N filter | `-n`, `--nions` | Show only the N highest-yield candidates (reference ion always included) |
| Highlight ions | (GUI field) | Space- or comma-separated ion names to highlight — user-selected only, never automatically assigned |

## Supported formats

- Spectra: `.npz` (configurable key mapping).
- Candidate lists: LISE++ output (`.lpp`).
- Export: `.ods` (candidate table), `simulation_result.out` (fixed-width
  text table).
- `.root` files are explicitly **not** supported — this is a deliberate
  restriction, not a bug; convert to `.npz` first.

## Troubleshooting

- **No network access on first run:** the AME2020 mass table is
  downloaded to `~/.ame/` on first use if not already cached there. If you
  have no network access, obtain `~/.ame/ame.data` from a machine that
  does, and copy the `~/.ame/` directory over.
- **Qt import crashes / binding-detection errors:** if you have both
  PyQt5 and PySide6 installed, some Qt-related tooling (notably
  `pytest-qt`, if you're running the test suite) can crash trying to
  auto-detect which binding to use. Ensure only one Qt binding is
  installed in your environment, or see `tests/conftest.py` for how this
  project works around it during testing.
- **GUI window doesn't appear / crashes with no display:** RionID needs a
  display server. On a headless machine, set `QT_QPA_PLATFORM=offscreen`
  for testing/scripting purposes (this disables the actual visible window
  — it's for automated checks, not for normal interactive use).

## Citation

If you use RionID, please cite it as described in
the repository's
[`CITATION.cff`](https://github.com/GSI-Nuclear-Astrophysics/RionID/blob/master/CITATION.cff).

## Limitations

- The polynomial correction is setting-specific: coefficients fitted for
  one ring configuration, optics, or cooling condition are not valid for
  another without re-validation.
- A quadratic correction order is justified only by the anchor residuals
  observed over the fitted interval — it is not a first-principles model,
  and extrapolation outside that interval is not supported.
- Candidate-list completeness bounds every assignment: an unlisted species
  cannot be identified.
- Unresolved blends and weak lines may remain ambiguous even with
  multi-harmonic, multi-detector comparison.
- RionID supports identification preceding a mass or lifetime measurement;
  it does not replace dedicated mass-calibration or lifetime-analysis
  procedures.
- This release performs no automatic or autonomous species assignment —
  see "Scope and non-goals" above.

## Acknowledgements

*   **Dr. RuiJiu Chen** for providing the C++ Time-of-Flight simulation code
    that inspired the backbone of this software.
*   **Dr. Shahab Sanjari** for guidance on software architecture and
    Schottky analysis.

## License

This project is licensed under the GNU General Public License v3.0. See
the [LICENSE](https://github.com/GSI-Nuclear-Astrophysics/RionID/blob/master/LICENSE)
file for details.
