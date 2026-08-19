# Open Scientific/Numerical Questions

Status: **not fixed, not physics-locked-in** — each item here is a
suspected pre-existing issue found during audit or tooling (not
introduced by any change in this repository's `joss-epja-prep`/JOSS-
packaging work), deliberately left as-is per this project's rule not to
silently alter numerical/physics-adjacent behaviour. Each entry states
the evidence, the potential impact, and a proposed validation — resolving
any of these requires your explicit decision.

## 1. `ImportData.calc_ref_rev_frequency`'s `gamma` is unbound if none of `brho`/`ke`/`gam` is provided

**Location:** `src/rionid/core.py`, `calc_ref_rev_frequency` (static method).
**Evidence:** `pyright --pythonversion 3.9 src/` reports `"gamma" is
possibly unbound (reportPossiblyUnboundVariable)`. The method's `if
brho: gamma = ...\nelif ke: gamma = ...\nelif gam: gamma = gam` has no
final `else`, unlike its caller `reference_frequency`, which does guard
the equivalent three-way check with `else: sys.exit(...)`.
**Potential impact:** `calc_ref_rev_frequency` is only reachable today via
`reference_frequency`, which already validates that at least one of
`fref`/`brho`/`ke`/`gam` is set before calling it — so this is not
reachable through the current call graph. It would become a live
`UnboundLocalError` risk only if `calc_ref_rev_frequency` were ever called
directly (it's a `@staticmethod`, technically public) without that
upstream guard.
**Proposed validation:** either add the same `else: raise ValueError(...)`
guard `reference_frequency` already has, or explicitly document that this
static method assumes a validated caller. No behaviour change to any
currently-working input.

## 2. `LISEreader._read`'s `file_start` is unbound if the `[Calculations]` marker is missing

**Location:** `src/rionid/external/lisereader/reader.py`, `_read`.
**Evidence:** `pyright` reports `"file_start" is possibly unbound` at two
use sites. `file_start` is set only inside `if '[Calculations]' in line:
file_start = i+1` within a `for` loop over the file's lines, then used
immediately after the loop with no fallback.
**Potential impact:** a malformed or unexpected LISE++ output file (missing
the `[Calculations]` section header) would raise `NameError` rather than a
clear, actionable error message about the file format.
**Proposed validation:** add an explicit check after the search loop
(`if file_start is None: raise ValueError(...)`) with a message naming the
expected `[Calculations]` marker. This is error-message quality, not a
numerical change — safe to fix once you approve it, but not done here.

## 3. `LISEreader.get_index`'s `element` is unbound if `self.data` is empty

**Location:** `src/rionid/external/lisereader/reader.py`, `get_index`.
**Evidence:** `pyright` reports `"element" is possibly unbound` (×2) in the
`else` branch reached when no matching entries are found — `element` is the
loop variable from the preceding `for i, element in enumerate(self.data):`,
unbound if `self.data` is empty.
**Potential impact:** an empty parsed LISE++ file would raise `NameError`
instead of the intended `ValueError` with a formatting hint.
**Proposed validation:** same class of fix as #2 — an explicit empty-data
guard with a clear message. Not done here.

## 4. `save_simulation_results`'s `harmonic`/`fre` rely on a loop-variable leak

**Location:** `src/rionid/gui/controller.py`, `save_simulation_results`.
**Evidence:** `pyright` reports both `"harmonic" is possibly unbound` and
`"fre" is possibly unbound`. The function's header-writing loop
(`for harmonic in harmonics: ...`) leaves `harmonic` bound via Python's
loop-variable-leaks-into-enclosing-scope behaviour, then the *last* value
of `harmonic` from that loop is used inside a *separate* loop further down
when `mode == 'Brho'`: `fre = mydata.srrf[i] * mydata.ref_frequency*harmonic`.
**Potential impact:** two things worth your attention, not just an unbound
risk: (a) if `harmonics` is ever empty, this is a real `NameError`; (b) even
when `harmonics` is non-empty, using only the *last* harmonic in the list
to compute every row's `fre` in Brho mode looks intentional-but-unobvious
for a multi-harmonic run — this is the kind of thing that's either a
correct simplification (Brho mode may only ever be called with one
harmonic in practice) or a latent output bug for multi-harmonic Brho runs,
and distinguishing those requires domain knowledge this audit doesn't have.
**Proposed validation:** confirm with the physics whether multi-harmonic
`-hrm` combined with Brho mode (`-b`) is a supported combination at all;
if not, an explicit early validation error would be clearer than the
current implicit last-value behaviour. If it is supported, the per-row
`fre` calculation likely needs the row's own harmonic, not the loop's last
one — a semantic question for you, not a bug this document is claiming to
have proven.

## 5. `rionid.__main__.run_controller`'s CLI entry point requires `-psim`/`--filep` despite argparse marking it optional

**Location:** `src/rionid/__main__.py`, `run_controller`.
**Evidence:** `run_controller` unconditionally calls
`ImportData._set_particles_to_simulate_from_file(particles_to_simulate)`
regardless of whether `-psim`/`--filep` was supplied, and `argparse`'s own
definition of `-psim` does not set `required=True`. Reproduced directly
(found while writing `docs/REPRODUCIBILITY.md`, Task 10 of the JOSS-
packaging plan): running `python3 -m rionid <any .npz> -r <ion> -ap <a>
-f <freq> -hrm <h>` without `-psim` raises `TypeError: expected str,
bytes or os.PathLike object, not NoneType` at
`src/rionid/external/lisereader/reader.py`'s `open(filename,
encoding="latin1")` — `LISEreader(None)` is constructed and immediately
tries to open `None` as a path. No `.lpp` fixture exists anywhere in this
repository, so this also cannot be reproduced end-to-end with a *working*
candidate file from a clean checkout — only worked around via direct
internal-API use (see `docs/REPRODUCIBILITY.md` §3).
**Potential impact:** every CLI invocation that omits `-psim` — which
argparse's own `--help` output implies is optional — crashes immediately
with an unhelpful `TypeError` instead of either (a) actually working
without a candidate list, or (b) a clear `argparse`-level error naming
`-psim` as required. This affects the CLI entry point's basic usability,
not any computed physics quantity: no code path here is capable of
returning a *wrong* numeric result, because it doesn't run at all.
**Proposed validation:** either mark `-psim`/`--filep` `required=True` in
the CLI argument parser (an honest CLI-contract fix, no numerical change)
or make `run_controller`/`_set_particles_to_simulate_from_file` handle
`None` gracefully (e.g., skip candidate-list loading entirely, simulating
only the reference ion) if that was ever meant to be a supported mode —
which of the two is correct depends on a product decision this document
isn't making.

## 6. `rionid.__main__`'s CLI entry point never supplies a ring `circumference`, so every reference-frequency mode crashes

**Location:** `src/rionid/__main__.py` (`main`/`run_controller`) and
`src/rionid/core.py` (`ImportData.__init__`, `_simulated_data`,
`calculate_brho_relativistic`).
**Evidence:** `ImportData.__init__` defaults `circumference=None` and
builds `self.ring = Ring("ESR", None)` when unset. `src/rionid/__main__.py`
has no `--circumference` CLI flag at all, and `run_controller` calls
`ImportData(ref_ion, alphap, filename=data_file)` without ever passing
one. `_simulated_data`'s non-`"brho"`-mode branch — the branch every CLI
run of `-f`/`-ke`/`-gam` takes, and even `-b` given the case mismatch
below — unconditionally calls `self.calculate_brho_relativistic(ref_moq,
ref_frequency, self.ring.circumference, harmonic)`, whose body computes
`v = actual_frequency * circumference`: `float * None` raises `TypeError:
unsupported operand type(s) for *: 'float' and 'NoneType'`. Reproduced
directly via an isolated call and via the full CLI path (see
`docs/REPRODUCIBILITY.md` §3's documented exception text; found alongside
item 5 above while writing that document). Separately, `_simulated_data`'s
branch check is `if mode == "brho":` (lowercase), while `__main__.py`/
`gui/inputs.py` pass capitalized mode strings (`"Frequency"`, `"Brho"`,
etc.) — so even explicit Brho-mode CLI runs take the same crashing
branch.
**Potential impact:** every CLI-driven simulation run crashes before
producing output, regardless of reference-frequency mode or arguments —
a usability/entry-point defect, not a computed-value error (nothing here
silently returns a wrong m/q, frequency, or yield; the run simply cannot
complete). The GUI path (`gui/controller.py::import_controller`) does
accept a `circumference` parameter from its own caller and is not
necessarily affected the same way — this item is specific to the
`python3 -m rionid` CLI entry point in `__main__.py`.
**Proposed validation:** add a `--circumference` CLI flag (mirroring the
GUI's own circumference input) and thread it through to `ImportData`,
plus reconcile the `mode` string casing between `__main__.py`'s callers
and `_simulated_data`/`save_simulation_results`'s checks (either
normalize case on read, or standardize the string values used
throughout). Both are CLI-usability/wiring fixes with no numerical
ambiguity — but out of scope for the JOSS-packaging documentation work
that discovered them.
