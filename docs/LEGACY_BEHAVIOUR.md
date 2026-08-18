# Legacy Behaviour (Phase 0 Audit)

Status: **audit only — no source changes**. This document records what the
codebase does *today*, as of commit `56ac07a` on branch `joss-epja-prep`
(pyproject version `8.0.0`), before any refactor, PID removal, or speed work
begins. It is the reference point for verifying that later changes are
numerically equivalent. See `docs/PUBLICATION_TRACEABILITY.md` for the
mapping to the EPJ A manuscript, `docs/AUTOMATIC_PID_REMOVAL_MAP.md` for the
automatic-PID inventory, and `docs/PERFORMANCE_BASELINE.md` for timings.

All line numbers refer to the current (pre-refactor) files under `src/rionid/`.

## Pipeline overview

```
file --> [1 ingest] --> [2 baseline] --> [3 normalise] --> [4 peak-detect]
                                                                  |
candidate list --> [5 mass/charge] --> [6 reference freq] --> [7 correction]
                                                                  |
                                              [8 harmonic projection] --> [10 top-N filter]
                                                                  |
                                                          [11 render] --> [12 export]
```
Stage 9 (automatic peak/species matching, "Quick PID") is documented
separately in `docs/AUTOMATIC_PID_REMOVAL_MAP.md` because it is out of
scope for this release.

## Stage 1 — Data ingestion

`ImportData.__init__` (`core.py:61`) dispatches on file extension via
`_get_experimental_data` (`core.py:178`) to `io.py`:

| Extension | Handler | Notes |
|---|---|---|
| `.csv` | `read_psdata` (`io.py:143`) | pipe-delimited, `dbm=False` reads column 1 |
| `.bin_fre/.bin_time/.bin_amp` | `handle_read_tdsm_bin` (`io.py:54`) | memory-mapped, time-averaged |
| `.npz` | `handle_spectrumnpz_data` (`io.py:103`) | keys default `arr_0`/`arr_1`, overridable via `io_params` (GUI: `KeySelectionDialog`, `gui/dialogs.py`) |
| `.root` | rejected | `core.py:194` raises `ValueError` — ROOT support was intentionally dropped (commit `04ae86e`) |

Loaded data is cached to `<basename>_cache.npz` (`_save_experimental_data`,
`core.py:297`) and reused unless `reload_data=True`
(`_load_experimental_data`, `core.py:303`).

**Dead/unreachable paths** (not automatic-PID, general cleanup candidates —
see `AUTOMATIC_PID_REMOVAL_MAP.md` decision items): `handle_tiqnpz_data`
(`io.py:72`) is imported into `core.py` but never called — the branch that
would call it is commented out (`core.py:188-192`). `handle_prerionidnpz_data`
(`io.py:124`) is defined but not imported/called anywhere.

## Stage 2 — Baseline removal (optional)

`NONPARAMS_EST.pls('BrPLS', ...)` (`baseline.py:23`) → `PLS.BrPLS`
(`baseline.py:65`), a Bayesian reweighted penalised-least-squares baseline
estimator (Wang et al., *Nucl. Sci. Tech.* 33:148, 2022). Applied twice in
`core.py` under `remove_baseline`: once in `_get_experimental_data`
(`core.py:196-204`) and again in the "NEW DATA PROCESSING BLOCK" in
`__init__` (`core.py:114-122`) — **the second application is the one that
actually reaches `self.experimental_data`** in the normal constructor flow;
the first (inside `_get_experimental_data`) only takes effect on the
`reload_data=True` path before the second block re-derives `freq, amp` from
`self.experimental_data`, so baseline subtraction is not applied twice in
practice, but the duplication itself is worth flattening during the refactor
(behaviour-preserving simplification, not a physics change).
Parameters: `psd_baseline_removed_l` (smoothness λ, default `1e6`),
fixed `ratio=1e-6`.

## Stage 3 — Log-safety clipping and normalisation

`core.py:124-136`: amplitude is floored at `1e-9` (`np.maximum`, comment says
"1e-9" but the code literal is `1e-29` at `core.py:127` — **discrepancy
between comment and code**, logged for `docs/OPEN_SCIENTIFIC_QUESTIONS.md`,
not fixed here) then divided by its max so the highest peak is `1.0`. This
means **displayed/exported amplitudes are relative, not absolute**, in every
current code path — an important unit/provenance fact for
`SCIENTIFIC_METHOD.md`.

## Stage 4 — Peak detection

`detect_peaks_and_widths` (`core.py:206-236`): `scipy.signal.find_peaks` with
`height = max(amp) * peak_threshold_pct`, `distance = min_distance`,
`prominence = height*0.2`, `width=1`. Peaks are then filtered to
`[matching_freq_min, matching_freq_max]` if set. Populates
`self.peak_freqs`/`self.peak_heights`; `self.peak_widths_freq` is always
zeroed (`core.py:236` — never actually computed from `find_peaks`' width
output, despite `peak_widths` being imported at `core.py:7` and unused).

## Stage 5 — Candidate definition and mass/charge

- `_parse_ref_ion` (`core.py:142-165`) accepts both `AAEl+QQ` and `AAElQQ+`
  reference-ion strings via regex, normalising to `AAElQQ+`.
- `_set_particles_to_simulate_from_file` (`core.py:313-318`) wraps
  `LISEreader` (`external/lisereader/reader.py`) to parse a LISE++ output
  file into `(element, aa, zz, nn, [charge_states], yield)` tuples, using the
  bundled AME2020 table (`external/barion/amedata.py`) to cross-reference
  each fragment.
- `_calculate_moqs` (`core.py:320-342`) has two call conventions: with an
  explicit `particles` list of `Particle` objects, or (the path actually used
  by the CLI/GUI) via `self.particles_to_simulate`, doing a **linear scan of
  the full AME table per candidate** (`core.py:334-336`) — see
  `PERFORMANCE_BASELINE.md` for the measured O(N×M) cost.
- Ionic mass (`Particle.get_ionic_mass_in_u`, `external/barion/particle.py:226`)
  subtracts the removed electrons' rest mass and adds the corresponding
  electron-binding-energy difference, via `AMEData.get_elbien`
  (`external/barion/amedata.py:197`, a `ZZ × charge-state` lookup table) —
  this is the code path implementing the manuscript's "atomic → ionic mass"
  correction (Sec. 2.1, `main.tex:198-200`). **Physics-locked; do not alter.**

## Stage 6 — Reference frequency

`reference_frequency` (`core.py:406-412`) supports four mutually exclusive
input modes — `fref` (Hz, direct), `brho` (Tm), `ke` (MeV/u), `gam`
(Lorentz γ) — routed through `calc_ref_rev_frequency` (`core.py:414-421`) and
the relativistic helpers `gamma_brho`/`gamma_ke`/`beta`/`velocity`
(`core.py:423-430`). `calculate_brho_relativistic` (`core.py:397-404`) is the
inverse direction (frequency → Brho), used when the GUI mode is `Frequency`
to report the corresponding Brho.

## Stage 7 — First-order model and polynomial correction (canonical)

`_calculate_srrf` (`core.py:344-356`) is **the single implementation** of the
manuscript's Eq. (1)/(2) forward model and Eq. (9)/(10) empirical residual
correction:

```python
self.srrf = array([1 - alphap*(moq[name]-moq[ref])/moq[ref] for name in moq])   # Eq. (2)
if correct:
    correction = polyval(array(correct), self.srrf * self.ref_frequency)         # Eq. (9): x = f^(0)
    self.srrf = self.srrf + correction / self.ref_frequency                      # Eq. (10)/ref_frequency
```
`correct = [A, B, C]` in `numpy.polyval` order (quadratic, linear, constant),
matching `main.tex:327-331` exactly. **Verified numerically in this audit**
against the manuscript's own worked example (harmonic-214 → revolution-space
→ harmonic-127 transform, Eqs. 15-21): reproduces every published digit, and
the code's "apply correction in revolution-frequency space, then multiply by
harmonic" order matches the direct harmonic-127 polynomial to `0.000e+00 Hz`
(script: see audit notes; recommend porting into a committed regression test
in Phase 1). **This function is the canonical correction and the single
place a future refactor's `correction.py` must reproduce exactly.**

`-c/--correct` (CLI) and the "Correction (a0\*x\*\*2 + a1\*x + a2)" field
(GUI, `gui/inputs.py:246-248`) are the only entry points; there is no
duplicate correction implementation elsewhere in the codebase (confirmed by
`polyval`/`correct` grep across `src/`).

## Stage 8 — Harmonic projection

`_simulated_data` (`core.py:358-395`) builds
`simulated_data_dict[str(harmonic)] = stack([freq, yield, name])` for each
requested harmonic, applying `Eq. (3)` (`F = h·f`) to the already-corrected
`srrf`. Yield is looked up by **a linear scan of `particles_to_simulate` per
candidate** (`core.py:373-381`) — the measured O(N²) hot path, see
`PERFORMANCE_BASELINE.md`.

## Stage 10 — Top-N display filter

`display_nions` (present independently in both `__main__.py:176-190` and
`gui/controller.py:124-163` — **two copies of the same logic**, a
simplification target). Sorts by yield descending, always re-includes the
reference ion even if outside the top N. Purely deterministic and
user-parameterised (`-n/--nions`); not automatic species selection.

## Stage 11 — Rendering

`CreatePyGUI` (`gui/plot.py`) is a `pyqtgraph` `QMainWindow`.
`plot_all_data` (`gui/plot.py:105-120`) always clears and fully redraws both
experimental and simulated layers — no incremental artist updates. Each
non-highlighted, non-reference ion label is a separate `pg.TextItem`
(`gui/plot.py:218-221`); the source code itself flags this as slow at scale
(`gui/plot.py:209`). Standard (non-highlight) lines for a harmonic are drawn
as one batched `PlotCurveItem` with `connect='pairs'`
(`gui/plot.py:239-255`) — already a reasonable vectorisation; the per-ion
`TextItem` labels are not. Log-Y axis is enabled by default
(`gui/plot.py:70`); amplitude floor `1e-9` is used again here for
log-safety, independent of the `1e-29` floor in `core.py` (two different
floor constants for the same purpose — flagged, not changed).

## Stage 12 — Export

- `write_arrays_to_ods` (`io.py:169`) — Name/freq/yield ODS export, triggered
  by `-w/--ods` or nothing in the GUI (no ODS control wired into
  `gui/inputs.py` — CLI-only today).
- `save_simulation_results` (`gui/controller.py:165-211`) — writes
  `simulation_result.out`, a fixed-width text table (ion, freq[Hz],
  yield[pps], m/q[u], mass[eV]) for **every** simulated ion, unfiltered by
  the top-N display filter (i.e. the export can contain more ions than are
  shown, if `nions` was used only for display).
- `save_matched_result` (`core.py:288-295`) — writes automatically-matched
  ions; see `AUTOMATIC_PID_REMOVAL_MAP.md` (remove).

## Known invariants to preserve across any refactor

1. Correction coefficient order `(A,B,C)` = (quadratic, linear, constant),
   applied as `Δf(x)=Ax²+Bx+C`, `x=f^{(0)}` (Hz), correction applied in
   revolution-frequency space *before* harmonic multiplication.
2. `alphap > 1` is silently reinterpreted as `γ_t` and converted via
   `alphap = 1/γ_t²` (`__main__.py:59-60`, `gui/controller.py:83`) — a
   real, intentional input-convenience convention, not a bug.
3. Reference ion is always force-included in top-N display/highlighting.
4. Amplitudes displayed/exported are baseline-and-max-normalised, not
   absolute, in every current path.
5. `.root` files are explicitly unsupported (raises, does not silently
   degrade).

## Legacy quirks noted but *not* fixed (physics/behaviour lock)

Per the audit's instructions, suspected issues are recorded here and in
`docs/OPEN_SCIENTIFIC_QUESTIONS.md` (to be created before Phase 2 edits),
not corrected silently:

- `core.py:127` floors amplitude at `1e-29` while the adjacent comment says
  "1e-9"; `gui/plot.py` separately floors at `1e-9`.
- `ImportData.__init__` always constructs `Ring('ESR', circumference)`
  (`core.py:77`) regardless of the actual facility, and never uses
  `Ring.gamma_t`/`Ring.get_alpha_p()` (`alphap` is taken directly as a
  constructor argument instead) — `Ring.get_ring_dict()`'s per-facility
  presets (ESR/CR/CSRe/CRYRING/CSR/HESR/TSR/RI-RING,
  `external/barion/ring.py:29-57`) are therefore unused by RionID.
- `AMEData.check_for_database` (`external/barion/amedata.py:83-96`) silently
  downloads the AME2020 table from `www-nds.iaea.org` on first use if no
  cache exists and no UI callback is supplied ("if no ui interface, just do
  it, don't ask") — a network side effect worth documenting in
  `REPRODUCIBILITY.md`.
- `save_simulation_results` always writes to the CWD-relative
  `simulation_result.out` (default arg), not the CLI's `-o/--outdir`.

## Dependency note (found during this audit, confirmed by grep)

`pyproject.toml` declares `requests = "^2.31.0"  # Required for Barion to
download mass tables`, but **`requests` is imported nowhere in `src/`** —
`external/barion/amedata.py:13` uses the standard-library `urllib.request`
for its AME-table download, not `requests`. This is a genuinely unused
runtime dependency. Flagged for removal in the same pass that extracts the
used `barion` routines into `masses.py` (see `REFACTORING_PLAN.md`), since
it directly touches the same file.

## Environment note (not a repository issue)

A stray `rionid==1.0.0` package is installed in the system Anaconda
`site-packages`, unrelated to this repository's `8.0.0` source tree. All
audit measurements in this document set used `sys.path` pointed explicitly
at `src/`, not the installed package. This will need `pip install -e .` in a
clean environment for Phase 5 verification, per the plan.
