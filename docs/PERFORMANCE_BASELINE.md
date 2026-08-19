# Performance Baseline (Phase 0 Audit)

Status: measurements taken against the **current, unmodified** code
(`src/rionid/core.py`, branch `joss-epja-prep`) via its public
`ImportData` API, on synthetic, non-confidential data. No source files were
changed to take these measurements.

## Hardware / environment

- Host: Linux 6.16.1-arch1-1, x86_64, 16 CPUs (`os.cpu_count()`).
- Python 3.9.7 (repo declares support for 3.9-3.12; measured on the low end
  of that range).
- No `Xvfb`/`xvfb-run` installed in this sandbox — worked around using
  PyQt5's built-in `QT_QPA_PLATFORM=offscreen` platform plugin instead,
  which runs the real Qt/`pyqtgraph` paint pipeline without any X server
  (confirmed working; see "GUI rendering" below). A few interaction-latency
  items still need a driven event loop and are listed under "Not measured."
- AME/NUBASE mass tables pre-cached at `~/.ame/` (no network access used).

## Method

`time.perf_counter()`, median of 3-5 repeats per call, each call invoking
only existing, unmodified public/semi-public methods of `ImportData`. Data:
a synthetic 200,000-point spectrum (`np.linspace` over a 1 MHz band + a few
Gaussian lines, seeded RNG) — clearly synthetic, not derived from any real
experiment, safe to redistribute as a fixture later. Candidate ion lists of
size N=10/100/2000 were built from **real** AME table rows (real Z/N/A/name,
fully-stripped, synthetic yield=1.0) rather than fabricated nuclides, to
keep the mass/charge arithmetic physically meaningful.

## Results

| Stage | N=10 | N=100 | N=2000 | Scaling |
|---|---:|---:|---:|---|
| Cold `import rionid.core` (once per process) | — | 254 ms | — | fixed |
| `AMEData()` table parse (cached, no network) | 168 ms | (same) | (same) | fixed per construction — **not currently cached across repeated `ImportData` runs** |
| `ImportData()` — load 200k-pt spectrum, baseline off, peak-detect | 130 ms | 130 ms | 131 ms | flat — dominated by npz I/O + `find_peaks`, independent of candidate count (candidates aren't loaded at construction) |
| `_calculate_moqs()` | 0.7 ms | 7.6 ms | 214 ms | **~O(N × 3558)** — linear AME-table scan per candidate (`core.py:334-336`) |
| `_calculate_srrf()`, no correction | 0.003 ms | 0.018 ms | 0.36 ms | O(N), cheap |
| `_calculate_srrf()`, with (A,B,C) correction | 0.013 ms | 0.029 ms | 0.38 ms | O(N), cheap — polynomial correction itself is not a bottleneck at any tested scale |
| `_simulated_data()`, 1 harmonic | 0.04 ms | 1.25 ms | **458 ms** | **~O(N²)** — nested per-candidate yield lookup (`core.py:373-381`) |
| `_simulated_data()`, 7 harmonics (123-129) | 0.12 ms | 1.61 ms | 466 ms | Same hot path; harmonic count adds only ~2% at N=2000, confirming the yield lookup (not the per-harmonic loop) dominates |
| `detect_peaks_and_widths()` re-run | 1.3-1.4 ms | ″ | ″ | flat, independent of candidate count (operates on the spectrum, not candidates) |

Raw script and console output are preserved in the audit trail (not
committed as a fixture yet — see Phase 1 note below).

## Reading the numbers

- The **polynomial correction itself (`_calculate_srrf`) is fast at every
  tested scale** (≤0.4 ms even at N=2000) — it is not a target for
  optimisation and must not be touched for performance reasons.
- Three genuine hot spots exist, **all pre-existing in the current code,
  all safe to fix without changing any output or export value**:
  1. `_calculate_moqs`'s per-candidate linear scan of the 3558-row AME table
     (`core.py:334-336`) — replaceable by a `{(name, aa): row}` dict built
     once per `AMEData` load, turning O(N×M) into O(N).
  2. `_simulated_data`'s per-`moq`-key linear scan of `particles_to_simulate`
     for yield lookup (`core.py:373-381`) — replaceable by a
     `{ion_name: yield}` dict built once, turning O(N²) into O(N).
  3. `plot_all_data`'s per-ion `pg.TextItem` label creation
     (`gui/plot.py:209,218-221`), which dominates wall-clock at N=2000
     (≈0.7-1.2 s, measured below) — a *display* optimisation (reuse/update
     artists, or an explicit, documented, exact level-of-detail rule at very
     high N), not a numerical one: `simulated_data_dict`/export values must
     stay byte-identical regardless of what's drawn, per the "never
     downsample/truncate to make display faster" rule.
  The first two are "efficient numeric/vectorised operations... [with]
  numerical results and error behaviour remain[ing] unchanged" per the speed
  rules; the third is a rendering change subject to the same output-identity
  requirement. All three are candidates for the next sub-project, pending
  your approval, with before/after values re-measured against this table.
- A fourth, non-algorithmic win: `AMEData()` re-parses the on-disk mass table
  (168 ms) on **every** `ImportData` construction, i.e. on every GUI "Run"
  click and (in the pre-removal code) on every element of a Quick-PID
  α_p×frequency grid that calls `_set_particles_to_simulate_from_file`
  again. Caching the parsed table for the process lifetime (immutable input,
  no invalidation complexity) removes this fixed cost from every run after
  the first.
- `ImportData()` construction (~130 ms, flat) is dominated by spectrum I/O
  and `scipy.find_peaks` on 200k points, not by anything candidate-related;
  this is the right order of magnitude to treat as "first usable spectrum
  render" prep cost baseline for a 200k-point input.

## GUI rendering (update: measured via `QT_QPA_PLATFORM=offscreen`)

Xvfb is not installed in this sandbox, but PyQt5's built-in **offscreen
platform plugin** (`QT_QPA_PLATFORM=offscreen`) constructs real
`QApplication`/`QMainWindow`/`pyqtgraph` objects and runs the actual paint
pipeline without any X server — confirmed working here, and usable in CI
without installing anything extra. This closes most of the "not measured"
gap from the first pass of this document.

| Stage | N=10 | N=100 | N=2000 |
|---|---:|---:|---:|
| `CreatePyGUI()` cold construction (once) | 34 ms | (same) | (same) |
| `plot_all_data()` — full clear + redraw (experimental + simulated) | 10 ms | 42 ms | **694-1216 ms** (high variance across repeats) |
| `reset_view()` (pan/zoom-equivalent view reset) | 0.15 ms | 0.15 ms | 0.15 ms — flat, independent of N |
| Peak RSS (cumulative across this single process, N=10→100→2000 in sequence) | 159 MB | 164 MB | 215 MB |

Reading: `plot_all_data` is cheap through N≈100 but becomes clearly
user-perceptible at N=2000 (roughly 0.7-1.2 s per redraw), consistent with
the source code's own comment that per-ion `pg.TextItem` creation
(`gui/plot.py:209,218-221`) "is still slow" at scale — this is a **separate,
additional** hot path from the two `core.py` hot spots above (moq lookup,
yield lookup), stacking on top of them. A single "Run" click at N=2000 on
this hardware costs roughly moq(214 ms) + simulate(458 ms) +
paint(700-1200 ms) ≈ **1.4-1.9 s end to end** — a real, perceptible lag,
and the clearest quantitative justification for the GUI speed work.
`reset_view()` (the cheapest available proxy for pan/zoom here, since no
real mouse-drag event loop was driven) shows no N-dependence, which is
expected — it only changes axis ranges, it doesn't redraw simulated items.

RSS grows monotonically across the N=10→100→2000 sequence run in the same
process (+86 MB total by N=2000), but this pass reused one growing `model`
object with an increasing candidate count each step, so growth is expected
from legitimately larger working-set data, not demonstrated to be a leak. A
proper leak test — repeated redraws **at fixed N**, watching whether RSS
keeps climbing rather than plateauing — is still open (see below).

## Not measured in this pass

- **True pan/zoom/tooltip latency under a driven mouse/wheel event loop**
  (this pass only called `reset_view()` as a proxy, and `mouse_moved`/
  cursor-label updates were not exercised) — offscreen Qt can drive
  synthetic `QTest` mouse events for this; left for the GUI-test sub-project
  rather than this baseline pass.
- **Leak-specific memory test** (repeated redraw at fixed N) — see above.
- **"Largest practical number of ions"** — 2000 was used as a stand-in for
  "large candidate list from a LISE++ fragmentation run"; if real use cases
  routinely exceed this, re-run with the actual scale before trusting the
  extrapolation blindly.

## What this baseline is for

Any Phase 2 change to `_calculate_moqs`, `_simulated_data`, or the AME-table
loading path must be re-measured with this exact script/methodology and
must show **identical output values** (same `moq`, `srrf`,
`simulated_data_dict` contents) at every N — only wall-clock time may
change. If a difference in output is detected, that is evidence to report
and stop on, per the task brief, not something to accept as a new default.

## Post-Wave-2a re-measurement

Status: measurements taken 2026-08-19 against the **current, optimized**
code on the same branch (`joss-epja-prep`), after all of Wave 2a's four
speed fixes landed. Same hardware (Linux 6.16.1-arch1-1, x86_64, 16 CPUs),
same Python (3.9.7), same methodology as Phase 0 above (`time.perf_counter()`
median of 5 repeats — the Phase 0 pass used 3-5; this pass used a
flat 5 throughout), same synthetic 200,000-point spectrum and real-AME-row
candidate lists at N=10/100/2000 (`tests.fixtures.synthetic_spectrum`).
Every timing run also asserted the *same output-shape invariants* the
regression suite checks (`len(moq) == len(yield_data) == N`, every
synthetic yield `== 1.0`, correct simulated-label count) — i.e. this was
not a speed-only script, output identity was checked inline as well as via
the dedicated tests below. The throwaway timing scripts used for this pass
were not committed (per the task brief); the numbers below are pasted
directly from their console output.

### Before/after: the four fixes

| # | Fix | Where | N=10 before → after | N=100 before → after | N=2000 before → after |
|---|---|---|---:|---:|---:|
| 1 | `_calculate_moqs` O(N×M) linear AME scan → O(N) dict lookup (`AMEData.lookup`/`._index`) | `masses.py:63-66`, `core.py:224-246` | 0.7 ms → 0.011 ms | 7.6 ms → 0.111 ms | **214 ms → 2.146 ms (~100×)** |
| 2 | `_simulated_data` O(N²) nested yield scan → O(N) `yield_by_name` dict | `core.py:262-302`, esp. 283-288 | 0.04 ms → 0.032 ms | 1.25 ms → 0.115 ms | **458 ms → 1.948 ms (~235×)** |
| 3 | `AMEData` table parse: re-parsed on every `ImportData`/`LISEreader` construction → process-lifetime cache shared by both call sites | `masses.py:726-741` (`get_ame_data()`); consumed by `core.py:219` and `external/lisereader/reader.py:15` | 168 ms/construction → **174.8 ms once, then 0.0002 ms/call** | (same, fixed cost) | (same, fixed cost) |
| 4 | Per-label `QFont("Arial", ...)` construction hoisted out of the per-ion redraw loop | `gui/plot.py:188` (built once before the loop, was previously built per label) | 10 ms → 10.1 ms | 42 ms → 41.6 ms | 694-1216 ms → 702.7-1719.3 ms (no material change — see note below) |

Fix 3's "before" number is actually a **conservative** restatement of Phase
0: `external/lisereader/reader.py`'s `LISEreader` was independently
constructing its own second `AMEData()` on every run (see the in-code
`NOTE` at `external/lisereader/reader.py:7-14`), a duplicate parse Phase 0
did not isolate — the real old per-run AME cost was closer to **2× the
measured 168 ms (~337 ms)** than the single 168 ms this document originally
recorded. Both `core.py` and `LISEreader` now share the one
`masses.get_ame_data()` cache, so the fixed cost is paid once per process
regardless of how many times either code path runs.

Fix 4's numbers show **no measurable improvement at any N** — this is
expected and consistent with the in-code performance note at
`gui/plot.py:162-178`: profiling during Task 9 found `plot_all_data`'s
N=2000 cost is dominated by PyQtGraph's `addItem`/`removeItem` scene-graph
reparenting overhead (~0.7-1.2 s of the ~1.4 s total), not by
`pg.TextItem`/`QFont` construction (~0.14 ms/item) — the QFont hoist
removes a real but small constant-factor cost per item, not the dominant
one. **The deferred TextItem-reuse optimization (diffing the previous vs.
new label set and reusing objects across redraws instead of
destroying/recreating them) remains open** — it requires visual QA this
sandbox cannot perform (no display server) and was deliberately not
attempted, per that same in-code comment and Task 9's write-up. This
re-measurement confirms that decision was correctly scoped: `plot_all_data`
at N=2000 is, within run-to-run variance, exactly as slow after Wave 2a as
before it.

### Full re-measured tables

`_calculate_moqs()` / `_simulated_data()` / `get_ame_data()` (headless, no Qt):

| Stage | N=10 | N=100 | N=2000 |
|---|---:|---:|---:|
| `get_ame_data()` first call (cold parse) | 174.8 ms | (same, once per process) | (same) |
| `get_ame_data()` subsequent calls (cache hit) | 0.0002 ms | 0.0002 ms | 0.0002 ms |
| `ImportData()` construction (200k-pt spectrum, peak-detect) | 133.7 ms | 133.3 ms | 133.1 ms — flat, as expected (unaffected by these fixes) |
| `_calculate_moqs()` | 0.011 ms | 0.111 ms | 2.146 ms |
| `_simulated_data()`, 1 harmonic | 0.032 ms | 0.115 ms | 1.948 ms |
| `_simulated_data()`, 7 harmonics | 0.112 ms | 0.467 ms | 8.070 ms |

`plot_all_data()` (offscreen Qt, `QT_QPA_PLATFORM=offscreen`), single harmonic per candidate:

| Stage | N=10 | N=100 | N=2000 |
|---|---:|---:|---:|
| `CreatePyGUI()` cold construction (once) | 33.8 ms | (same) | (same) |
| `plot_all_data()` — full clear + redraw | 10.1 ms | 41.6 ms | 702.7-1719.3 ms, median 1403.2 ms (high variance, unchanged from Phase 0) |
| `reset_view()` | 3.1 ms | 6.7 ms | 84.5 ms |
| Peak RSS (cumulative, N=10→100→2000 in one process) | 161 MB | 165 MB | 230 MB |

One honest discrepancy versus Phase 0, recorded rather than smoothed over
per the task brief: Phase 0 reported `reset_view()` as flat (~0.15 ms)
"independent of N," but this pass measured a clear N-dependence (3.1 ms →
6.7 ms → 84.5 ms). `reset_view()`'s own code
(`gui/plot.py:320-328`) is unchanged by Wave 2a — it only calls
`setXRange`/`setYRange` — so this is most plausibly the cost of Qt/
PyQtGraph re-clipping/repainting the now-larger scene (up to 2000
`TextItem`s + curves) triggered by the range change, not a regression this
plan introduced, and not one of the four fixes' targets. Flagged here for
visibility rather than silently reconciled with the older number.

### Output-identity confirmation

Every timing run in this pass carried inline assertions
(`len(model.moq) == N`, `len(model.yield_data) == N`, every yield `== 1.0`,
and — for the GUI pass — exactly `N` rendered labels), all of which held at
every N. These are the same invariants the committed regression suite
checks unconditionally:

- `tests/test_analysis.py::test_calculate_moqs_output[10|100]` and
  `test_simulated_data_yield_lookup_output[10|100]` /
  `test_simulated_data_yield_lookup_output_n2000` (the `-m slow` test) —
  golden output shape/value checks for fixes #1 and #2 above.
- `tests/test_masses.py::test_ionic_mass_and_moq[...]` and
  `test_get_ame_data_is_cached` — golden ionic-mass/m-q values (unaffected
  by the O(1) lookup rewrite) and the cache-identity check (`get_ame_data()`
  returns the *same* object on repeated calls) for fix #3.
- `tests/test_gui_smoke.py::test_plot_simulated_data_label_count_and_text`
  — label count/text identity for fix #4's `QFont` hoist.

All of the above pass on this branch as of this re-measurement (see the
full-suite run recorded in the Task 10 report) — every fix is confirmed
**faster with identical output**, exactly as Phase 0 required.
