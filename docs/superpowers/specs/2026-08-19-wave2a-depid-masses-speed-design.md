# Wave 2a: Regression tests, automatic-PID removal, `masses.py` extraction, baseline-removal removal, and safe speed fixes

Status: approved design, pre-implementation. Written after the Phase 0 audit
(`docs/LEGACY_BEHAVIOUR.md`, `docs/PUBLICATION_TRACEABILITY.md`,
`docs/AUTOMATIC_PID_REMOVAL_MAP.md`, `docs/PERFORMANCE_BASELINE.md`,
`REFACTORING_PLAN.md`), which this spec assumes as read and does not repeat.

## Why "Wave 2a" and not the full `REFACTORING_PLAN.md` restructuring

`REFACTORING_PLAN.md` describes a target 10-module layout
(`correction.py`/`calibration.py`/`analysis.py`/`species.py`/`models.py`/
`configuration.py`/`cache.py`/`exceptions.py`, plus GUI file renames). Doing
that simultaneously with deletion, extraction, and performance work would
produce one large, hard-to-review, hard-to-bisect change, and — with no
local display server for automated GUI verification — only one opportunity
to manually smoke-test the result instead of several smaller ones. This
sub-project ("Wave 2a") instead does the substantive, behaviour-affecting
work first, keeping `core.py`/`gui/inputs.py`/`gui/app.py`/`io.py` at their
current names, just smaller. The full module split ("Wave 2b") is deferred
to its own later sub-project, where it is a pure move with no logic change
— safest done last, verified by the test suite this wave builds.

## Scope of this wave

1. Build the public regression-test suite, against **current, unmodified**
   code, as characterisation tests.
2. Remove the automatic-PID surface (`docs/AUTOMATIC_PID_REMOVAL_MAP.md`
   "Remove" list, in full).
3. Extract the used subset of vendored `external/barion` into a new
   `src/rionid/masses.py`; delete `external/barion/`; drop the unused
   `requests` dependency from `pyproject.toml`.
4. Remove the baseline-subtraction feature (`baseline.py`, BrPLS) and all
   its call sites — **added to scope in this design's approval round**, see
   "Baseline removal" below.
5. Apply the three speed fixes identified in `PERFORMANCE_BASELINE.md`
   (`_calculate_moqs` dict-index, `_simulated_data` yield dict-index,
   `AMEData` process-lifetime cache) plus the `plot_all_data` label-creation
   investigation.
6. Re-measure performance against `PERFORMANCE_BASELINE.md`; re-run the full
   suite; confirm zero automatic-PID / zero baseline-removal surface
   remains in UI/config/docs.

`external/lisereader/` is untouched. No equation, unit, coefficient
ordering, or the canonical `_calculate_srrf` correction changes in any way.

## Baseline removal (added in design review)

You asked to also remove baseline subtraction and its dependency, and
flagged that ion labels are sometimes mispositioned when a baseline has
been removed. Scope:

- Delete `baseline.py` (`NONPARAMS_EST`, `PLS.BrPLS`) entirely.
- Delete both call sites in `core.py`: the "NEW DATA PROCESSING BLOCK" in
  `ImportData.__init__` (`core.py:114-122`) and the duplicate call inside
  `_get_experimental_data` (`core.py:196-204`) — `LEGACY_BEHAVIOUR.md`
  already notes the second one is redundant with the first in the normal
  constructor flow, so removing the feature also removes that duplication
  for free.
- Delete `remove_baseline`/`psd_baseline_removed_l` from `ImportData.__init__`'s
  signature, `import_controller`'s signature (`gui/controller.py`), and the
  GUI (`remove_baseline_checkbox`, `psd_baseline_removed_l_edit`, and their
  `load_parameters`/`save_parameters` keys in `gui/inputs.py`).
- **No `pyproject.toml` dependency actually needs removing**: `baseline.py`
  only imports from `scipy` (`scipy.special`, `scipy.sparse`,
  `scipy.sparse.linalg`), and `scipy` stays required regardless
  (`scipy.signal.find_peaks` in `core.py`'s peak detection is unaffected and
  stays). Flagging this now so "and its dependency" isn't silently
  misread as a `pyproject.toml` change that doesn't actually apply.
- README.md's "Built-in baseline subtraction (BrPLS) and peak detection"
  feature bullet and the `--remove_baseline`/`--peak_threshold_pct` CLI-flag
  documentation lines are removed. **Confirmed during this design pass**:
  these two flags are documented in `README.md`'s "Arguments" section but
  **do not exist in `__main__.py`'s argparse at all**, and `run_controller`
  never threads `remove_baseline`/`peak_threshold_pct` through to
  `ImportData` on the CLI path either — baseline removal was already
  GUI-only and unreachable from the CLI. This independently confirms the
  removal is self-contained (no CLI surface to touch) and separately fixes
  a pre-existing README/CLI documentation mismatch.
- The reported mispositioned-label symptom is not being root-caused —
  removing the feature makes the question moot rather than something to
  debug. Recorded here for provenance only, not carried into
  `docs/OPEN_SCIENTIFIC_QUESTIONS.md` (that file is for physics/statistics
  issues in retained functionality, and this one is being removed, not
  fixed-in-place).
- This is a product/UI simplification decision made directly by you (the
  domain owner), not a "suspected physics issue" being silently patched —
  recorded here as the authorization for the change.

## Test suite (built first; must pass against current code before any deletion begins)

All tests live under `tests/`, use `pytest`, and run in public CI (no real
experimental data, no restricted paths — confirmed not needed this round).

- **`tests/test_correction.py`** — the manuscript's harmonic-214→127
  numeric example (`main.tex` Eqs. 15-21), ported from the Phase 0
  verification script into an assertion-based golden test:
  `_calculate_srrf`'s correction step must reproduce every published digit
  of the worked example, and the "correct in revolution space, then
  multiply by harmonic" path must equal the direct harmonic-127 polynomial
  to machine precision (`rtol=1e-12` — pure `float64` arithmetic, no
  measurement uncertainty involved, so this is tight by design).
- **`tests/test_masses.py`** — a handful of named, real nuclides (e.g.
  ¹²C⁶⁺, ⁷²Ge³²⁺, ⁷⁴As³³⁺, ⁷⁶Se³⁴⁺ — chosen because they appear in the
  manuscript's own worked example) with ionic mass / m-q computed through
  `external.barion.particle.Particle` (present until step 3 deletes it).
  Once `masses.py` exists, the same assertions are re-pointed at it and
  must match to `rtol=1e-12`; **`external/barion/` is only deleted once this
  test is green against `masses.py`**, per `REFACTORING_PLAN.md`'s risk #2.
- **`tests/test_analysis.py`** — `_calculate_moqs`/`_simulated_data` at
  small N (fast, always run) and a `@pytest.mark.slow` N=2000 variant
  (opt-in, matches `PERFORMANCE_BASELINE.md`'s scale), asserting **identical
  `moq`/`simulated_data_dict` contents** before and after the dict-index
  speed fixes — this test is written and run against the current O(N²)
  implementation first, establishing the golden output, then re-run
  unmodified after the fix lands.
- **`tests/fixtures/`** — the synthetic-spectrum generator and
  real-AME-row candidate builder from the Phase 0 profiling scripts,
  committed as reusable, clearly-synthetic fixtures (docstring states
  plainly this is not experimental data).
- **`tests/test_gui_smoke.py`** (`pytest-qt`, `QT_QPA_PLATFORM=offscreen`)
  — construct `MainWindow`; assert Quick-PID widgets/attributes
  (`setup_quick_pid`, `quick_pid_script`, `alphap_min_edit`, etc.) and
  baseline-removal widgets (`remove_baseline_checkbox`,
  `psd_baseline_removed_l_edit`) are **absent** after removal; assert
  `run_script` still emits `visualization_signal` and the plot updates;
  assert `overlay_sim_signal`/`_stop_quick_pid`/`onPlotClicked` no longer
  exist. This is state/interaction-based (per the brief's instruction to
  avoid screenshot-only GUI tests), not pixel comparison.
- Floating-point comparisons use explicit `rtol`/`atol` chosen per quantity:
  `1e-12` relative for pure-arithmetic identities (correction formula,
  mass/m-q), exact equality (`==`) for anything that should be bit-identical
  dict/array output before vs. after a *pure* refactor (dict-index fixes),
  never a blanket default tolerance.

`pytest-qt` is added to `[tool.poetry.group.dev.dependencies]` only (no
runtime/shipped dependency change).

## Speed fixes (applied after the suite is green and PID/baseline removal has landed)

1. `_calculate_moqs`: build `{(name, aa): row}` once from the AME table
   (in `masses.py`, alongside the table loader) instead of scanning per
   candidate. O(N×3558) → O(N).
2. `_simulated_data`: build `{ion_name: yield}` once from
   `particles_to_simulate` instead of scanning per `moq` key. O(N²) → O(N).
3. Cache the parsed AME/NUBASE table for the process lifetime (module-level,
   immutable input, no invalidation logic needed) so repeated `ImportData`
   construction within one GUI session doesn't re-parse it every time.
4. `plot_all_data`'s per-ion `pg.TextItem` creation: investigate a
   same-output speed-up (e.g. reuse/update existing `TextItem`s across
   redraws instead of destroy-and-recreate, since positions/labels for
   unchanged ions don't need new Qt objects). **If no safe, exact-output
   equivalent is found, this is reported as not done, not forced** — per
   the "never alter display to fake speed" rule, a level-of-detail approach
   would need your explicit sign-off first and is not assumed here.

Every fix is re-measured against `docs/PERFORMANCE_BASELINE.md`'s exact
methodology (same synthetic fixture, same N=10/100/2000 points), with output
identity confirmed by `tests/test_analysis.py` before the timing comparison
is trusted.

## Error handling

Bare `except Exception: print(...)`, `sys.exit(...)`, and print-then-raise
patterns are converted to a small `exceptions.py` hierarchy (e.g.
`RionIDError`, `UnsupportedFileFormatError`, `ReferenceIonParseError`) **only
in functions already being edited** for PID removal, baseline removal, or
the speed fixes — this wave does not do a blanket exception-handling sweep
across untouched code (that belongs to Wave 2b or later, where those files
are being restructured anyway).

## Acceptance criteria

- [ ] Every test in "Test suite" above green.
- [ ] Zero automatic-PID surface remains (re-check against every "Remove"
      row in `AUTOMATIC_PID_REMOVAL_MAP.md`).
- [ ] Zero baseline-removal surface remains in code, GUI, README, or
      `parameters_cache.toml` schema.
- [ ] `external/barion/` deleted; `masses.py` exists and
      `tests/test_masses.py` passes against it; `requests` removed from
      `pyproject.toml`.
- [ ] `PERFORMANCE_BASELINE.md` re-measured, showing improvement at
      N=10/100/2000 with **identical** `moq`/`simulated_data_dict` output
      (or a documented reason a fix wasn't safe to apply).
- [ ] `ruff check`/`ruff format --check` clean on every touched file.
- [ ] You manually smoke-test the GUI before merge (no display in this
      sandbox for a final human-eyes check).

## Explicitly out of scope for this wave

- The full `correction.py`/`calibration.py`/`analysis.py`/`species.py`/
  `models.py`/`configuration.py`/`cache.py`/`exceptions.py` module split
  and GUI file renames (Wave 2b, separate sub-project).
- Uncertainty propagation for the polynomial correction (Eq. 22 gap) —
  product-scope decision, not a refactor.
- The Eq. (8) harmonic-overlap-count calculation.
- Any change to `_calculate_srrf`'s arithmetic, coefficient ordering, or
  units.
- The DOI mismatch (`docs/PUBLICATION_TRACEABILITY.md`) — flagged only, per
  your earlier answer.
- `bibliography.bib`/`thebibliography` reconciliation — Phase 4 (manuscript)
  work.
