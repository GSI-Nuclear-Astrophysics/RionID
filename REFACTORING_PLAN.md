# Refactoring Plan

Status: **plan only — describes the target architecture and migration
sequence; no code has been moved yet.** Written after the Phase 0 audit
(`docs/LEGACY_BEHAVIOUR.md`, `docs/PUBLICATION_TRACEABILITY.md`,
`docs/AUTOMATIC_PID_REMOVAL_MAP.md`, `docs/PERFORMANCE_BASELINE.md`), which
this plan assumes as read. Execution of any step below is its own
sub-project and gets its own brainstorm → design → approval cycle before
code changes begin, per the agreed sequencing.

## Design goals (from the brief, restated as constraints)

1. Physics/numerics separated from GUI, I/O, plotting, configuration.
2. One canonical polynomial-correction implementation (already true today —
   `ImportData._calculate_srrf` is the only place `polyval`+`correct` is
   used — the refactor's job is to keep it that way while relocating it).
3. Zero automatic-PID surface (see `AUTOMATIC_PID_REMOVAL_MAP.md`).
4. GUI never blocks on repeated expensive fits/scans/disk reads (Quick PID's
   `QApplication.processEvents()`-driven blocking loop is the one place this
   happens today, and it is being deleted, not fixed, per the removal map).
5. Every parameter explicit, typed, validated, documented.
6. No empty abstractions added just to match a template tree.

## Target layout

```
src/rionid/
├── app.py            # NEW: thin composition/entry point (QApplication + MainWindow wiring)
├── gui/
│   ├── main_window.py    # renamed from app.py's MainWindow class (composition root moves up to app.py)
│   ├── controller.py      # bridges GUI inputs -> analysis/correction/calibration modules; stays
│   ├── panels.py          # renamed from inputs.py: widget construction + parameter (de)serialization only
│   ├── dialogs.py         # stays (KeySelectionDialog)
│   └── spectrum_view.py   # renamed from plot.py: PyQtGraph widget (stays GUI-owned; see rationale below)
├── io.py              # stays: file-format readers/writers (read_psdata, handle_*_data, write_arrays_to_ods)
├── models.py          # NEW: typed data objects (SpectrumData, CandidateIon, CorrectionCoefficients,
│                       #      SimulationResult) replacing today's loose dict/tuple bags
├── configuration.py   # NEW: parameters_cache.toml load/save, extracted out of gui/panels.py
├── correction.py      # NEW: the canonical polynomial correction, extracted verbatim from
│                       #      ImportData._calculate_srrf's correction step
├── calibration.py     # NEW: reference-frequency determination (reference_frequency,
│                       #      calc_ref_rev_frequency, calculate_brho_relativistic, gamma/beta/velocity helpers)
├── analysis.py         # NEW: pure scientific calculations (mass/charge orchestration via _calculate_moqs,
│                       #      harmonic projection via _simulated_data, detect_peaks_and_widths)
├── species.py          # NEW: reference-ion parsing, highlight-ion parsing, LISE++ candidate ingestion
├── plotting.py? — see rationale below (NOT created; spectrum_view.py stays under gui/)
├── export.py           # NEW: save_simulation_results, write_arrays_to_ods re-export, ODS/table export
├── cache.py            # NEW: spectrum npz cache (existing) + AMEData/candidate-list cache (new,
│                       #      addresses the PERFORMANCE_BASELINE.md re-parse finding)
├── masses.py           # NEW: extracted-and-trimmed from external/barion (ionic mass/m-q incl.
│                       #      electron-binding correction, AME/NUBASE table loader, unit constants,
│                       #      ring-circumference holder) — see "barion extraction" below
├── exceptions.py       # NEW: RionIDError hierarchy replacing bare except/sys.exit/print-and-raise
├── __main__.py         # stays: CLI entry point, now calling analysis/calibration/correction directly
└── external/
    └── lisereader/      # unchanged: vendored LISE++ parser (GPL-3.0). external/barion/ is REMOVED
                        # as a vendored dependency — see "barion extraction" below (decision made
                        # 2026-08-18, AUTOMATIC_PID_REMOVAL_MAP.md decisions #1-2)
```

### Deliberate deviations from the example tree (and why)

- **No top-level `plotting.py`.** The template assumes a plotting module
  that can produce figures independent of a GUI event loop. RionID's
  `CreatePyGUI` is inherently a `QMainWindow` — its rendering methods
  (`plot_all_data`, `plot_experimental_data`, `plot_simulated_data`) are not
  meaningfully separable from Qt without inventing a parallel headless
  plotting path nobody asked for. Keeping it under `gui/spectrum_view.py`
  avoids an empty abstraction. If a future need arises for headless figure
  export (e.g. for the paper's Fig. 1 regeneration), that would be a new,
  explicitly-requested `export.py` function that reuses the same data
  objects, not a `plotting.py` split of the widget.
- **`app.py` becomes the entry point**, with the current `gui/app.py`
  `MainWindow` class renamed into `gui/main_window.py` and `app.py`
  reduced to `main()` + `QApplication` setup — this satisfies "app.py: GUI
  composition/entry point only" literally while keeping `MainWindow` itself
  inside `gui/`.

## Legacy → new mapping

| Current | New | Notes |
|---|---|---|
| `core.py: ImportData.__init__` (data load + baseline + normalise + peak-detect) | `io.py` (load) + `analysis.py` (baseline/normalise/peak-detect) + `models.py: SpectrumData` | Constructor becomes an orchestration function, not a 140-line `__init__`. |
| `core.py: _parse_ref_ion`, `_parse_highlight_ions` | `species.py` | Pure parsing, easily unit-tested in isolation. |
| `core.py: _calculate_moqs` | `analysis.py` + AME-lookup dict from `cache.py` | Same output; O(N×M) → O(N) per `PERFORMANCE_BASELINE.md`. |
| `core.py: _calculate_srrf` | `correction.py` (correction step) + `calibration.py` (reference-frequency step) | Split along the manuscript's own Eq. (2) vs Eq. (9)/(10) boundary — reference-frequency determination is calibration, the polynomial residual is correction. **The arithmetic itself does not change.** |
| `core.py: _simulated_data` | `analysis.py` + dict-based yield lookup from `models.py: CandidateIon` list | Same output; O(N²) → O(N). |
| `core.py: compute_matches`, `chi2`, `match_count`, `save_matched_result` | **Deleted** | Per `AUTOMATIC_PID_REMOVAL_MAP.md`. |
| `core.py: calculate_brho_relativistic`, `reference_frequency`, `calc_ref_rev_frequency`, `gamma_brho/gamma_ke/beta/velocity` | `calibration.py` | Moved verbatim. |
| `__main__.py: display_nions` + `gui/controller.py: display_nions` (duplicated) | `analysis.py: filter_top_n_ions` (single copy) | Behaviour-preserving de-duplication. |
| `gui/inputs.py` (573 lines: widgets + quick-PID + parameter cache) | `gui/panels.py` (widgets + parameter I/O calls into `configuration.py`) | Quick-PID code deleted, not moved. |
| `gui/inputs.py: CollapsibleGroupBox` (unused) | Deleted, pending your confirmation (removal map, decision #4) | |
| `gui/app.py: MainWindow` | `gui/main_window.py` | `overlay_sim_signal`/`overlay_simulation` deleted (removal map item 8). |
| `gui/plot.py` | `gui/spectrum_view.py` | Rendering logic unchanged; incremental-update work (Phase 2 speed pass) lands here later, as its own reviewed change. |
| `io.py` | `io.py` (unchanged location) | `handle_tiqnpz_data`/`handle_prerionidnpz_data` **deleted** (removal-map decision #3, resolved). |
| `baseline.py` | `analysis.py` (imported, not merged — BrPLS is a distinct, citable algorithm; keeping it in its own file is also acceptable and lower-diff) | **To decide during implementation**, not a physics question. |
| `external/barion/particle.py: Particle.__init__` (zz/nn match), `get_ionic_mass_in_u`, `get_ionic_moq_in_u` | `masses.py` | Extracted verbatim (arithmetic unchanged, including the electron-binding correction). |
| `external/barion/amedata.py: AMEData` table loader, `to_mev`/`to_u`/`to_kg`, `get_elbien`, `CC`/`UU`/`EE`/`ME` | `masses.py` | Extracted verbatim; this is the only part of `amedata.py` (837 lines) actually reachable from RionID (confirmed by call-graph in `LEGACY_BEHAVIOUR.md`). |
| `external/barion/ring.py: Ring` | `masses.py` (reduced to a circumference holder — `gamma_t`/`get_alpha_p()`/`get_ring_dict()` presets are unused by RionID, per `LEGACY_BEHAVIOUR.md`) | Extracted, trimmed. |
| `external/barion/particle.py`: `identify_range`/`identify_search_range`/`identify_region`/`identify_search_region`/`get_unknown_moq_from_freq`/`get_unknown_rev_freq_from_moq`/`get_nuclides_freqs`/`get_all_in_all`/`get_isotopes`/`get_isotones`/`get_isobars` | **Deleted, not extracted** | Removal-map decision #1/#2, resolved — you are also an owner of `barion` upstream, so this is a scoping decision, not a vendoring risk. |
| `external/barion/` (whole subpackage) | **Removed** | `pyproject.toml` loses the vendored-`barion` framing entirely; only `masses.py`'s extracted subset remains, as RionID's own code. |
| `external/lisereader/` | Unchanged | Out of scope for this decision (LISE++ parsing only; not part of the barion question). Flag if this reading is wrong. |

## Migration approach

**Characterisation-tests-first**, since none exist today:

1. Write the Phase 1 regression suite (starting with the correction-formula
   and harmonic-invariant tests already specified by the manuscript and
   verified manually in this audit) against the **current, unrefactored**
   code. Every golden value is cross-checked against the manuscript's own
   worked numbers where one exists (Eqs. 9-21), not blindly snapshotted from
   current output — this avoids freezing an undiscovered bug as the new
   "correct" answer.
2. Only once that suite is green does any move/rename/split begin.
3. Each module extraction in the mapping table above is a pure move (copy
   function body, update imports, delete original) with no logic change,
   verified by re-running the same suite after each step.
4. Automatic-PID deletion (removal map "Remove" list) happens as its own
   commit(s) after the suite shows the remaining parameter-driven workflow
   does not depend on any of it (confirmed already in this audit by call-
   graph inspection — `compute_matches`'s only caller is `quick_pid_script`).
5. Performance changes (AME-table caching, dict-based lookups) are applied
   last, each re-measured against `docs/PERFORMANCE_BASELINE.md` with
   identical-output verification before/after.

## Acceptance criteria

- [ ] One canonical polynomial-correction call site (`correction.py`), used
      by both CLI and GUI — no duplicate implementation anywhere.
- [ ] Zero automatic-PID surface in UI, config, docs, or tests (verified
      against every "Remove" row in `AUTOMATIC_PID_REMOVAL_MAP.md`).
- [ ] Phase 1 regression suite green before and after every migration step.
- [ ] `PERFORMANCE_BASELINE.md`'s three identified hot spots improved with
      byte-for-byte identical output at N=10/100/2000, or the change is
      rejected/reported per the numerical-equivalence rule.
- [ ] `ruff check`, `ruff format --check`, `pyright` clean on the new layout.
- [ ] GUI manually smoke-tested by you before merge (this sandbox has no
      display server; automated GUI-launch verification is not possible
      here — see `PERFORMANCE_BASELINE.md` "Not measured").

## Risks

1. **No pre-existing tests** means "characterisation test" is new work, not
   a snapshot of trusted behaviour — mitigated by anchoring golden values to
   the manuscript's own numbers wherever one exists, and flagging (not
   silently accepting) any golden value that has no independent source.
2. **`masses.py` extraction must be checked against real AME/NUBASE data**,
   not just re-read: the electron-binding-energy lookup (`get_elbien`) and
   the AME table parser (`fortranformat`-based, `amedata.py`) are the two
   places most likely to have a subtle transcription error if copied by
   hand. The Phase 1 regression suite should include a golden test that
   loads a few known nuclides through both the pre-extraction (`external.
   barion`, still present in the current commit) and post-extraction
   (`masses.py`) code paths and asserts identical mass/m-q output, before
   `external/barion/` is deleted.
3. **No local GUI display** — every GUI-facing change in this environment
   needs a manual verification pass by you before it's trusted; this plan
   does not claim GUI correctness from static reading alone.
4. **DOI/CITATION.cff mismatch and manuscript TODO placeholders**
   (`PUBLICATION_TRACEABILITY.md`) are outside this refactor's scope but
   block any "release-ready" claim in Phase 5 until resolved.

## Explicitly out of scope for this plan

- Adding uncertainty propagation for the polynomial correction (Eq. 22 gap,
  `PUBLICATION_TRACEABILITY.md`) — real product-scope work requiring a
  design decision, not a refactor.
- Adding the Eq. (8) harmonic-overlap-count calculation to the package
  (currently only in the paper) — same reason.
- Any change to `_calculate_srrf`'s arithmetic, coefficient ordering, or
  units.
