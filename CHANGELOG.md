# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this
project uses [Semantic Versioning](https://semver.org/).

## Unreleased

- Summarize unavailable LISE++ masses once per file and export the omitted
  records with their reasons, preserving the previous duplicate-record policy.
- Keep saved-run replay independent of the live AME table, including LISE++
  parsing, and preserve the original import diagnostics on replay.

- Add manual ion candidates and charge states independently of LISE++.
- Search the complete AME2020 mass-to-charge catalogue with configurable
  relative tolerance and charge bounds; retain candidate provenance and
  extrapolation flags.
- Filter and highlight isotopes or elements across charge states and harmonics
  without recalculating the model.
- Propagate mass-input uncertainties to frequency bands using the corrected
  model derivative or an explicit local Hz/100 keV conversion.
- Render large overlays by viewport with bulk markers and a bounded label
  count, while retaining complete candidate tables.
- Support the same candidate and uncertainty settings in the CLI and saved
  configurations. Allow overlays without a LISE++ file.

## [9.0.1] - 2026-08-21

### Fixed

- Removed the obsolete Python `<3.13` package-metadata restriction so
  RionID can be installed on Python 3.13 and 3.14.

## [9.0.0] - 2026-08-21

Breaking changes to the public API and CLI. Earlier history (4.0.0-8.0.0)
is not reconstructed here from git log; see the git tag history for that.

### Removed

- **Automatic ion identification ("Quick PID")**: the GUI panel that
  scanned a grid of (reference frequency × α_p) values and silently
  overwrote user-entered values with an automatic best-match result.
  Removing it keeps RionID's workflow deterministic and expert-guided.
- **Baseline subtraction (BrPLS)**: the "Remove Baseline" feature and its
  `remove_baseline`/`psd_baseline_removed_l` parameters, per a deliberate
  product-scope decision (not a physics fix).
- **Peak utilities**: peak detection, peak plotting, peak threshold and
  matching-frequency controls were removed from the supported workflow.
- Dead/unreachable code: `handle_tiqnpz_data`/`handle_prerionidnpz_data`
  I/O handlers, and several
  unused GUI widgets tied to the removed features above.

### Changed

- Two algorithmic hot paths are now O(N) instead of O(N×M)/O(N²):
  candidate mass-to-charge lookup and simulated-data yield lookup.
- The AME2020 mass table is now cached for the process lifetime
  instead of being re-parsed on every simulation run.
- `pyproject.toml` migrated to PEP 621 metadata; `pip install -e ".[dev]"`
  now works (it could not before).

### Added

- A public pytest regression suite (27 default tests plus one slow test) — there were none before
  this release cycle. Covers the polynomial-correction arithmetic, the
  extracted mass-table arithmetic, GUI state/wiring after the removals
  above, and the optimized analysis paths' output identity.
- `ruff`/`pyright`/`pre-commit`/CI configuration.

### Fixed

- A stale-mass-table double-parse: `external/lisereader`'s LISE++ reader
  used to construct its own independent copy of the AME table on top of
  the one `ImportData` already loads — both now share one cached load.
