# Performance Baseline (Phase 0 Audit)

Status: measurements taken against the **current, unmodified** code
(`src/rionid/core.py`, branch `joss-epja-prep`) via its public
`ImportData` API, on synthetic, non-confidential data. No source files were
changed to take these measurements.

## Hardware / environment

- Host: Linux 6.16.1-arch1-1, x86_64, 16 CPUs (`os.cpu_count()`).
- Python 3.9.7 (repo declares support for 3.9-3.12; measured on the low end
  of that range).
- No display server (`Xvfb`/`xvfb-run` not installed in this sandbox) — see
  "Not measured" below for what this excludes.
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
- Two genuine algorithmic hot spots exist, **both pre-existing in the
  current code, both safe to fix without changing any output value**:
  1. `_calculate_moqs`'s per-candidate linear scan of the 3558-row AME table
     (`core.py:334-336`) — replaceable by a `{(name, aa): row}` dict built
     once per `AMEData` load, turning O(N×M) into O(N).
  2. `_simulated_data`'s per-`moq`-key linear scan of `particles_to_simulate`
     for yield lookup (`core.py:373-381`) — replaceable by a
     `{ion_name: yield}` dict built once, turning O(N²) into O(N).
  Both are "efficient numeric/vectorised operations... [with] numerical
  results and error behaviour remain[ing] unchanged" per the speed rules —
  candidates for the next sub-project, pending your approval, with
  before/after values re-measured against this table.
- A third, non-algorithmic win: `AMEData()` re-parses the on-disk mass table
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

## Not measured in this pass (needs a display / needs the user's machine)

- **Cold GUI startup, first paint, pan/zoom/tooltip latency** — no `Xvfb` in
  this sandbox. The computational path above (data load → moq → correction
  → simulate) is display-independent and fully covered; only the Qt/
  `pyqtgraph` paint loop is not. Recommend re-running this section on a
  machine with a display (or after adding `Xvfb` to CI) before/after the
  Phase 2 GUI work, using the same `ImportData` calls feeding
  `CreatePyGUI.plot_all_data`.
- **Memory (RSS) during repeated navigation** — not instrumented this pass
  (would need `tracemalloc`/`psutil` wrapped around the same calls above,
  plus a real paint loop for the GUI-retention question). Add in Phase 1.
- **"Largest practical number of ions"** — 2000 was used as a stand-in for
  "large candidate list from a LISE++ fragmentation run"; if real use cases
  routinely exceed this, re-run with the actual scale before trusting the
  O(N)/O(N²) extrapolation blindly.

## What this baseline is for

Any Phase 2 change to `_calculate_moqs`, `_simulated_data`, or the AME-table
loading path must be re-measured with this exact script/methodology and
must show **identical output values** (same `moq`, `srrf`,
`simulated_data_dict` contents) at every N — only wall-clock time may
change. If a difference in output is detected, that is evidence to report
and stop on, per the task brief, not something to accept as a new default.
