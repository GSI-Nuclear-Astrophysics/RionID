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

The GUI can read legacy local settings from `parameters_cache.toml`.
That file is ignored by Git. A sanitized template is provided as
`parameters_cache.example.toml`; copy it locally to prepare starting settings.
Use the visible saved-run controls below to retain a configuration and results;
running a simulation no longer overwrites the legacy settings cache.

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

## Saving configurations, tables and replaying a run

The **Saved runs: configuration and tables** panel is visible at the top of the
GUI. Choose an **Output parent** folder and optionally enter a **Folder label**
(for example, `203Pt` or `setting_A`).

- **Run** computes the overlay in a background thread.
- **Save config + tables** saves the current settings and matching results.
  When the form differs from the last completed run, it saves settings only.
- **Save results bundle (config + tables)** saves the last completed run with
  the settings that actually produced it, even if the form has since changed.
- **Load config…** restores settings without computing.
- **Replay config…** restores settings, reruns the simulation, updates the
  overlay, and checks the saved input fingerprints and numerical checkpoints.

For example, saving on 6 October 2026 with label `203Pt` creates:

```text
2026-10-06_RionID_CONFIG_203Pt/
    config.json
    candidates.csv
    candidates.ods
    simulation_result.out
```

The text table name can be changed using **Text table name**; it always remains
inside this folder. All tables include the ion label, harmonic, radiofrequency,
revolution frequency, relative yield, ionic mass and mass-to-charge ratio. Rows
include all calculated candidates and flag which are displayed after the top-N
filter. Relative yields are normalized LISE++ yields, optionally scaled for
plotting; they are not measured ion counts or production rates.

Repeated saves of unchanged settings and results reuse the same folder in the
current GUI session. Changed runs or labels create a new folder, with `_02`,
`_03`, etc. added when necessary. Existing snapshots are not overwritten. Saved
folders inside a checkout are ignored by Git.

The JSON snapshot records all simulation settings, the NPZ array-key mapping,
input-file SHA-256 fingerprints, the software and NumPy versions, and numerical
checkpoints. Completed runs also retain the resolved ionic mass inputs used by
the frequency calculation. Replay uses those saved inputs rather than querying
a possibly changed mass table. Settings-only configurations use the current
mass inputs when run and have no stored numerical results to verify.

When **Reload Data Cache** is unchecked and an existing `*_cache.npz` is used,
the snapshot records that cache as the effective spectrum input, including its
actual array keys. Replay reads that exact spectrum rather than silently
switching to the original file. Saved-run operations do not create or overwrite
source caches and do not copy raw spectra or LISE++ files into the output folder.
They currently support NPZ spectra or existing NPZ caches, matching the GUI's
experimental-data workflow.

If a recorded input has moved, Replay asks for its new location and verifies the
same content. Changed inputs are rejected. Frequency checkpoints use an absolute
`1e-5 Hz` comparison tolerance; this is a numerical replay tolerance, not a
measurement uncertainty. Version differences and unmatched numerical results
are reported separately.

Replay is also available through Python:

```python
from rionid.session import replay_configuration, save_bundle

report = replay_configuration("2026-10-06_RionID_CONFIG_203Pt/config.json")
assert report.verified, (report.warnings, report.mismatches)
folder, files = save_bundle(
    "outputs", "replayed", report.completed.config, report.completed
)
```

## Additional candidates and display

Expand **Additional ions, AME search and display** to add ions such as
`205Au76+, 205Au75+, 201Ir76+`. The LISE++ file is optional when manual ions
or an AME search supply the candidates. Duplicate ions are merged, and the
reference is included automatically. These additions use AME2020 ground-state
masses; they do not supply isomer excitation energies.

When a LISE++ candidate has no AME2020 mass, RionID omits that record and
reports one summary for the file. The saved configuration retains the import
counts and reasons; `skipped_candidates.csv` lists every omitted record.
Duplicate transport records keep the first matching record, as in earlier
versions. Frozen replay restores the saved ionic masses and import diagnostics
without querying the live AME table.

To explore nearby mass-to-charge ratios, enter `205Au78+` under **Search AME
around m/q of** and set the half-width, for example `1` for ±1%. The search
considers every AME2020 nuclide and allowed positive charge state, optionally
restricted by the charge range. It includes AME extrapolations and labels them
in the exported tables. A candidate satisfying this interval is a possible
frequency overlay, not an identification or a production prediction.

**Show only** and **Ions to highlight** accept exact ions (`205Au78+`),
isotopes (`205Au`, all available charge states) and elements (`Au`). Separate
selectors with spaces or commas. They apply across every harmonic. Changing
these fields updates the cached overlay without recalculating the mass model.
A filter can reveal matching candidates beyond the top-N yield selection.

Manual and AME-search candidates without LISE++ yields use dash-dot markers
at a guide height; their table yield is empty. Hover over a label to inspect
its source and mass-input uncertainty. Candidate lines are drawn in batches.
Large overlays render the current frequency window and at most 300 labels;
zooming refreshes the labels while all candidates remain in the saved tables.

### Mass-input frequency uncertainty

**Show ±1σ mass-input frequency bands** displays the uncertainty caused by
AME2020 mass inputs. By default, the derivative of the corrected frequency
model propagates target and reference mass errors, with the shared reference
handled consistently. The bands are not peak widths or a complete calibration
uncertainty: they exclude frequency-reference, correction-model and
electronic-binding uncertainties. AME extrapolation errors remain estimates
from the evaluation.

An optional empirical conversion can replace that derivative. For example,
enter `50` Hz per 100 keV and harmonic `214` if that conversion was determined
at the 214th harmonic. RionID rescales it for each displayed harmonic. In an
older configuration that uses an RF reference and harmonic *ratios*, supply
the corresponding multiplier instead (usually `1` at the central RF). This
option propagates the target mass error through an approximate local scale;
it does not determine or fit that scale. Its own uncertainty is held fixed,
and the conversion should be used only in the region where it was calibrated.
Distinct nuclide mass errors are treated as independent because the AME table
does not supply their covariance; charge states of the same nuclide share their
mass error. Press **Run** after changing the scale.

The same workflow is available from the CLI:

```bash
rionid-cli spectrum.npz -psim candidates.lpp \
  -r 204Hg77+ -ap 0.512 -f 1932523.3644859814 -hrm 210 211 212 213 214 \
  --add-ions '205Au76+,205Au75+' \
  --search-ion 205Au78+ --moq-tolerance-percent 1 --charge-min 70 \
  --display-filter Au --highlight-ions 205Au --show-uncertainties \
  --frequency-scale-hz-per-100kev 50 --uncertainty-harmonic 214 \
  --frequency-key frequency --amplitude-key amplitude -s
```

Omit `-psim` for a manual/search-only overlay. Use `-w -o outputs` to export
an ODS candidate table including sources, extrapolation flags and propagated
frequency errors. New settings and resolved mass inputs are retained in saved
GUI configurations, while older saved configurations remain readable.

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

- Spectra: `.npz` (configurable key mapping), `.csv` (pipe-delimited),
  `.bin_fre`/`.bin_time`/`.bin_amp` (TDSM binary triples).
- Candidate lists: LISE++ output (`.lpp`).
- Export: `.ods` (candidate table), `simulation_result.out` (fixed-width
  text table).
- `.root` files are explicitly **not** supported — this is a deliberate
  restriction, not a bug; convert to `.npz`/`.csv` first.

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
[`CITATION.cff`](CITATION.cff).

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
the [LICENSE](LICENSE) file for details.
