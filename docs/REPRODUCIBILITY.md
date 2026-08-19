# Reproducibility

Exact, runnable commands to reproduce this repository's key numerical
claims from a clean checkout, using only public, synthetic, or already-
committed data — none of these need real experimental data.

## 1. Reproduce the polynomial-correction golden check

This reproduces the manuscript's harmonic-214 → revolution-space →
harmonic-127 coefficient-transform worked example (`RionID-EPJA/main.tex`
Eqs. 9-23) directly against the shipped code:

```bash
pip install -e ".[dev]"
pytest tests/test_correction.py -v
```
Expected: `1 passed`. The test's own docstring/comments walk through the
derivation; read `tests/test_correction.py` alongside
`docs/PUBLICATION_TRACEABILITY.md` for the full manuscript mapping.

## 2. Reproduce the extracted-mass-table golden check

Confirms `rionid.masses`'s ionic mass/mass-to-charge calculation for four
named nuclides (including the manuscript's own E143 example: ⁷²Ge³²⁺,
⁷⁴As³³⁺, ⁷⁶Se³⁴⁺):

```bash
pytest tests/test_masses.py -v
```
Expected: `5 passed`. Requires network access on first run only, to
populate `~/.ame/` with the AME2020 mass table (cached afterward).

## 3. Reproduce a full simulation run on synthetic data

No real spectrum file needed — this generates one from
`tests/fixtures/synthetic_spectrum.py` (documented there as clearly
synthetic, not derived from any experiment):

```bash
python3 -c "
import sys
sys.path.insert(0, 'src')
from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum
make_synthetic_spectrum('synthetic_demo.npz')
print('wrote synthetic_demo.npz')
"
```
Expected: prints `wrote synthetic_demo.npz` and creates the file. Verified.

**A note on the packaged `python3 -m rionid <datafile>` CLI here:** as
currently wired, `rionid.__main__.run_controller` unconditionally calls
`ImportData._set_particles_to_simulate_from_file`, which requires a real
LISE++ candidate-list file via `-psim`/`--filep`; no such file is public
or synthetic-committed in this repository, so a CLI invocation without
`-psim` fails immediately with
`TypeError: expected str, bytes or os.PathLike object, not NoneType`
(`src/rionid/external/lisereader/reader.py`, `open(filename, ...)`).
Independently, `rionid.__main__` has no `--circumference` flag and never
passes one to `ImportData`, so `ImportData.ring.circumference` stays
`None`; `_simulated_data`'s Brho calculation then multiplies a frequency
by `None`, raising
`TypeError: unsupported operand type(s) for *: 'float' and 'NoneType'`
for every reference-frequency mode (`-f`/`-b`/`-ke`/`-gam`) alike — this
reproduces independently of the `-psim` issue. Both were confirmed by
direct reproduction while writing this document. Neither is fixable here
without touching `src/rionid/`, which is out of scope for this
documentation task.

What *is* reproducible with only public/synthetic, already-committed
data is the same underlying simulation engine and result-writing code
the CLI calls, driven the way `tests/test_analysis.py` drives it: real
AME2020 rows (real Z/N/A/element name; only the candidate selection and
per-candidate yield are synthetic) via
`tests/fixtures/synthetic_spectrum.build_ame_candidates`, in place of a
LISE++ file, with an explicit ring circumference (ESR, 108.36 m) supplied
directly since the CLI does not expose one:

```bash
python3 -c "
import sys
sys.path.insert(0, 'src')
from tests.fixtures.synthetic_spectrum import build_ame_candidates
from rionid.core import ImportData
from rionid.masses import get_ame_data
from rionid.gui.controller import save_simulation_results
from numpy import argsort

ame = get_ame_data()
candidates = build_ame_candidates(ame.ame_table, 750)
assert ('Ge', 72) in [(c[0], c[1]) for c in candidates]

mydata = ImportData('72Ge+32', alphap=0.189, filename='synthetic_demo.npz',
                     reload_data=True, circumference=108.36)
mydata.ame = ame
mydata.ame_data = ame.ame_table
mydata.particles_to_simulate = candidates
mydata._calculate_moqs()
mydata._calculate_srrf(fref=1930000, correct=[2.55222702e-7, -0.985167644, 950690.215])
mydata._simulated_data(harmonics=[127.0], mode='Frequency')

sort_index = argsort(mydata.srrf)
save_simulation_results(mydata, 'Frequency', [127.0], sort_index, filename='simulation_result.out')
print(f'wrote simulation_result.out with {len(mydata.nuclei_names)} candidates')
"
```
Expected (verified, deterministic): prints
`wrote simulation_result.out with 750 candidates` and writes
`simulation_result.out` — a fixed-width table of candidate ion
name/frequency/yield/m-q/mass, e.g. the reference ion's own line:
```
72Ge32+        1929995.7047597999            9.6667e-01     2.247018202880 66978687329.326
```
(m/q = 2.247018202880 u matches `tests/test_masses.py`'s independently
verified value for ⁷²Ge³²⁺ to displayed precision; the frequency sits
close to, not exactly at, the nominal 1930000 Hz reference because the
polynomial correction from §1 has been applied on top of it, as
intended.)

## 4. Reproduce the performance baseline

The exact synthetic-data methodology behind `docs/PERFORMANCE_BASELINE.md`
(both the Phase 0 numbers and Wave 2a's post-fix re-measurement) uses
`tests/fixtures/synthetic_spectrum.make_synthetic_spectrum`/
`build_ame_candidates` with `rionid.masses.get_ame_data()`, timed via
`time.perf_counter()` medians at N=10/100/2000 candidates — see that
document's "Method" section for the full, literal methodology; the
`@pytest.mark.slow` test in `tests/test_analysis.py`
(`pytest tests/ -m slow -v`) exercises the N=2000 scale as a correctness
(not timing) check and is the closest committed artifact to that
methodology.

## 5. What is *not* reproducible from this repository alone

`RionID-EPJA/main.tex`'s Figures 1, 3, and 4 (the broadband E143 spectrum,
the dual-pickup contaminant comparison, and the calibration-residual plot)
require the real E143 experimental dataset, which is not distributed in
this repository and is subject to the E143 collaboration's data policy —
see the manuscript's own Data Availability section. Figure 2 (the
harmonic-overlap-count staircase) is not reproducible from the package at
all, since the software does not compute that quantity — see
`docs/SCIENTIFIC_METHOD.md` §4.
