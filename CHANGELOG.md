# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this
project uses [Semantic Versioning](https://semver.org/).

## [9.0.0] - Unreleased

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
  I/O handlers, the vendored `barion` library's unused automatic-
  identification methods (`identify_range` and related), and several
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
