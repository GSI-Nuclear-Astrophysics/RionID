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
