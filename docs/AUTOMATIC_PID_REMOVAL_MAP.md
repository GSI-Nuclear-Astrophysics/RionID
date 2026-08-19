# Automatic-PID Removal Map (Phase 0 Audit)

Status: **all "Remove" items below are confirmed removed as of Wave 2a**
(re-verified 2026-08-19, Task 10 — see the Status column in the "Remove"
table and `docs/PERFORMANCE_BASELINE.md`'s companion re-measurement pass).
This document originated as an audit-only inventory before any deletion
happened, per the original task brief: every PID-related
file/function/UI item/config/test/dependency/claim, classified as
**Remove**, **Retain (rationale)**, or **Decision needed**. That removal
work is now complete — `src/rionid/gui/inputs.py:setup_quick_pid`,
`quick_pid_script`, `src/rionid/core.py:compute_matches`,
`save_matched_result`, and every other item in the "Remove" table below no
longer exist in the shipped code. This is confirmed both by direct
`grep`/`hasattr` checks (see the "Confirming zero automatic-PID surface
remains" section at the end of this document) and by the dedicated
regression tests in `tests/test_gui_smoke.py`
(`test_quick_pid_surface_is_absent`, `test_compute_matches_is_removed`,
`test_overlay_quick_pid_wiring_is_removed`, `test_collapsible_group_box_is_removed`,
`test_baseline_removal_surface_is_absent`, `test_baseline_module_is_removed`),
all of which pass on the current branch.

## What "automatic PID" means here, concretely

The manuscript (`RionID-EPJA/main.tex`) states throughout that RionID
"separates calculation from scientific acceptance" and is explicitly **not**
"an autonomous classifier" (abstract, conclusions). The shipped code
currently violates that in exactly one self-contained feature, wired end to
end: the **"Quick PID"** GUI panel, which scans a grid of
(reference frequency × α_p) values, computes a χ² match between simulated
lines and detected experimental peaks for each grid point
(`ImportData.compute_matches`), and **silently overwrites** the user's
reference-frequency field, α_p field, and `highlight_ions` selection with
the automatically "best" result. This is the entire automatic-PID surface;
it does not touch the polynomial-correction physics (§ verified separately
in `PUBLICATION_TRACEABILITY.md`).

## Remove

| # | Item | Location (pre-removal) | Why | Status |
|---|---|---|---|---|
| 1 | `setup_quick_pid()` | `gui/inputs.py:261-297` | Builds the "Quick PID" `QGroupBox` (alpha range, freq range, "Run Quick PID" button). | ✅ Removed in Wave 2a |
| 2 | `quick_pid_script()` | `gui/inputs.py:422-537` | The scanning algorithm itself: nested loop over experimental peaks × α_p grid, calls `compute_matches`, silently overwrites `value_edit`/`alphap_edit`, emits the result as if user-confirmed. | ✅ Removed in Wave 2a |
| 3 | `alphap_min_edit`, `alphap_max_edit`, `alphap_step_edit`, `fref_min_edit`, `fref_max_edit`, `threshold_edit` widgets, and their `load_parameters`/`save_parameters` keys (`alphap_min`, `alphap_max`, `alphap_step`, `fref_min`, `fref_max`, `threshold`) | `gui/inputs.py:88-93,119-124,231-236,265-290` | Config surface that exists only to feed #2. | ✅ Removed in Wave 2a |
| 4 | `ImportData.compute_matches()` | `core.py:238-286` | Nearest-neighbour peak↔simulated-line χ² matcher; its only caller is #2. Also the function that overwrites `self.highlight_ions` with automatically matched species (`core.py:284`) — the clearest violation of "user provides desired ions; software calculates and displays them." | ✅ Removed in Wave 2a |
| 5 | `ImportData.chi2`, `ImportData.match_count` attributes | `core.py:94-95` | State written only by #4. | ✅ Removed in Wave 2a |
| 6 | `ImportData.save_matched_result()` | `core.py:288-295` | Exports the automatically-matched ion list. Confirmed **dead** even today (no caller anywhere in the GUI), but conceptually part of the automatic-matching feature — remove rather than wire up. | ✅ Removed in Wave 2a |
| 7 | `matched_result_edit` field | `gui/inputs.py:87,118,254,256` | Loaded/saved to `parameters_cache.toml` but never connected to any action — a "half-supported flag" for #6. | ✅ Removed in Wave 2a |
| 8 | `overlay_sim_signal` (`pyqtSignal`) + `MainWindow.overlay_simulation()` | `gui/inputs.py:30`, `gui/app.py:75,95-110` | Docstring: "used primarily during the 'Quick PID' scan to show visual feedback... without reloading the heavy experimental data every frame." Once #2 is gone, nothing emits this signal. The normal full-refresh path (`visualization_signal` → `MainWindow.update_visualization`) is unaffected and stays. | ✅ Removed in Wave 2a |
| 9 | `RionID_GUI._stop_quick_pid`, `onPlotClicked()`, and the `MainWindow` wiring `visualization_widget.plotClicked.connect(rion_input.onPlotClicked)` | `gui/inputs.py:37,336-338,487,494,500`, `gui/app.py:78` | Stop-flag for #2's loop, set by clicking the plot. **Do not** remove `CreatePyGUI.plotClicked` itself (`gui/plot.py:37,58`) — it also drives the legitimate manual "Pick" coordinate feature (see Retain list). Only this specific subscriber goes. | ✅ Removed in Wave 2a (`CreatePyGUI.plotClicked` itself retained, confirmed still present) |
| 10 | README.md / `docs/index.md` "Automated Matching: Includes Quick PID logic to scan α_p and Reference Frequency to find the best match (χ² minimization)" bullet | `README.md`, `docs/index.md` | Unsupported-capability claim for this release; contradicts the manuscript's own "not an autonomous classifier" statement. | ✅ Removed in Wave 2a |

Verified 2026-08-19 by direct `grep -rniE
"quick.?pid|compute_matches|overlay_sim_signal" src/ README.md
docs/index.md` (Task 10 of
`docs/superpowers/plans/2026-08-19-wave2a-depid-masses-speed.md`)
returning zero matches, and by `tests/test_gui_smoke.py`'s
`test_quick_pid_surface_is_absent`, `test_compute_matches_is_removed`, and
`test_overlay_quick_pid_wiring_is_removed`, all passing on this branch.

Everything in this list is self-contained: removing it deletes one GUI panel,
one algorithm, and its dedicated config/signal plumbing, and does not touch
`_calculate_srrf` (the canonical correction), `_calculate_moqs`,
`_simulated_data`, `detect_peaks_and_widths`, or any CLI argument.

## Retain (legitimate deterministic calculation or visual aid — not automatic PID)

| Item | Location | Rationale |
|---|---|---|
| `detect_peaks_and_widths()` (`scipy.find_peaks`) | `core.py:206-236` | Deterministic peak-finding *on the experimental spectrum itself*; does not match against simulation or assign species. Feeds the "Detected Peaks" markers and gives the user candidate values to manually pick from. |
| `matching_freq_min_edit`/`matching_freq_max_edit` + "Pick" buttons | `gui/inputs.py:212-229` | A deterministic, user-set frequency window for peak detection/display — not species assignment. **Distinct** from Quick PID's `fref_min`/`fref_max` (item 3 above), which only ever fed the automatic scan grid. |
| `enterPlotPickMode`/`_onPlotPicked`, `CreatePyGUI.plotClicked` | `gui/inputs.py:311-329`, `gui/plot.py:37,58` | Manual cursor-to-value picking — an explicitly permitted visual aid, not automatic assignment. Retained in full; only the Quick-PID-specific subscriber (item 9) is removed. |
| `-c/--correct` / `correction_edit` (canonical polynomial correction) | `core.py:344-356`, `gui/inputs.py:246-248` | Core physics, entirely user-parameterised. Untouched. |
| `-n/--nions`, `display_nions()` | `__main__.py:176-190`, `gui/controller.py:124-163` | Deterministic "top N by yield" display filter; **N and the ranking criterion (yield) are fixed and user-set**, not a species match/rank produced by comparing to experimental data. Matches the brief's "user-selected number or list of candidate ions" requirement directly. |
| `highlight_ions` field | `gui/inputs.py:173`, `core.py:73,167-171` | Retained as a pure user-input display selector. Only its automatic mutation inside `compute_matches` (item 4) is removed; after that, nothing ever writes to it except the user. |

## Decisions (resolved 2026-08-18)

| # | Item | Location | Resolution | Status |
|---|---|---|---|---|
| 1 | `Particle.identify_range`, `identify_search_range`, `identify_region`, `identify_search_region`, `get_unknown_moq_from_freq`, `get_unknown_rev_freq_from_moq`, `get_nuclides_freqs` | `external/barion/particle.py:392-474` | **Remove.** You are also an owner of the upstream `barion` library, so vendoring/re-vendoring risk does not apply the way it would for an arbitrary third party. Decision: **drop the `external/barion` vendored subpackage as a dependency entirely**; extract only the routines RionID actually calls (ionic-mass/m-q calculation including the electron-binding correction, AME/NUBASE mass-table loading and unit constants, and a minimal ring-circumference holder) into RionID's own source tree. `identify_*`/`get_unknown_*`/`get_nuclides_freqs` are not among the extracted routines and are dropped. | ✅ Removed in Wave 2a — `external/barion/` deleted entirely; only `external/lisereader/` remains, plus the extracted subset in `src/rionid/masses.py` |
| 2 | `Particle.get_all_in_all`, `get_isotopes`, `get_isotones`, `get_isobars`, `Ring.get_ring_dict()` | `external/barion/particle.py:83-171`, `external/barion/ring.py:29-57` | **Remove**, same resolution as #1 — not part of the extracted subset RionID uses. | ✅ Removed in Wave 2a (deleted with all of `external/barion/`) |
| 3 | `handle_tiqnpz_data`, `handle_prerionidnpz_data` | `io.py:72-101,124-141` | **Remove.** Confirmed dead code; no re-wiring requested. | ✅ Removed in Wave 2a |
| 4 | `CollapsibleGroupBox` | `gui/inputs.py:546-574` | **Remove**, unused, unrelated to PID. | ✅ Removed in Wave 2a |

None of these are automatic-PID in themselves (only #1 overlaps with automatic
identification, and only incidentally, as an inert unused algorithm) — they
are recorded here because their disposition was decided in the same review
pass as the PID removal. The extraction scope for #1/#2 is: `Particle.qq`,
`Particle.__init__` matching-by-(zz,nn), `get_ionic_mass_in_u`,
`get_ionic_moq_in_u`, and `AMEData`'s table loader, unit-conversion statics
(`to_mev`/`to_u`/`to_kg`), `get_elbien`, and the physical constants
(`CC`/`UU`/`EE`/`ME`) — i.e. exactly the call graph confirmed reachable from
`core.py` in `LEGACY_BEHAVIOUR.md`. **Scope note**: this decision covers
`external/barion` only; `external/lisereader` (the LISE++ parser) was not
part of the question and is unaffected — flag if that reading is wrong.

## Dependencies

No Python dependency in `pyproject.toml` exists solely for automatic PID —
`compute_matches`/`quick_pid_script` use only `numpy` (already required for
the core physics) and Qt signals already used elsewhere. **No dependency
removal is needed or possible from this change alone.** (Confirmed still
true post-removal — no `pyproject.toml` change was required by this work.)

## Tests

✅ Done. `tests/test_gui_smoke.py` (added during Wave 2a) is the regression
suite this section originally called "upcoming": it asserts the full
"Remove" list is absent by `hasattr`/`grep`-equivalent checks
(`test_quick_pid_surface_is_absent`, `test_compute_matches_is_removed`,
`test_collapsible_group_box_is_removed`, `test_overlay_quick_pid_wiring_is_removed`)
and that the retained signal path still works
(`test_visualization_signal_still_updates_the_plot`). All pass on this
branch as of the Task 10 re-verification (`python3 -m pytest tests/ -v`).

## Documentation/claims to update once removal lands

- ✅ `README.md` "Features" bullet (item 10 above) — updated, no
  `quick.?pid`-matching text remains (verified by grep, see below).
- ✅ `docs/index.md` (identical bullet) — updated identically.
- `RionID-EPJA/main.tex` already correctly claims *no* automatic
  classification — no manuscript text changes were needed for this reason
  specifically (out of scope for this document; any other manuscript
  review is tracked elsewhere).
- `parameters_cache.toml` (local, untracked) simply stopped being written
  the removed keys once the code was gone; no action was needed, confirmed.

## Confirming zero automatic-PID / baseline-removal surface remains

Re-run 2026-08-19 (Task 10 of
`docs/superpowers/plans/2026-08-19-wave2a-depid-masses-speed.md`):

```bash
grep -rniE "quick.?pid|compute_matches|overlay_sim_signal|remove_baseline|psd_baseline" src/ README.md docs/index.md
```

Result: **zero matches.** The only remaining occurrences of these terms in
the repository are in historical planning/audit documents (this file,
`docs/PERFORMANCE_BASELINE.md`, `docs/LEGACY_BEHAVIOUR.md`,
`docs/PUBLICATION_TRACEABILITY.md`, `REFACTORING_PLAN.md`,
`docs/superpowers/specs/` and `docs/superpowers/plans/`) and in
`tests/test_gui_smoke.py`, where the names appear only as strings being
asserted *absent* via `hasattr`/`importlib` checks — neither is shipped
surface. The baseline-subtraction feature (`remove_baseline_checkbox`,
`psd_baseline_removed_l_edit`, the `rionid.baseline` module, BrPLS) was
added to this document's original PID-focused scope during design review
(see `docs/superpowers/specs/2026-08-19-wave2a-depid-masses-speed-design.md`)
and is confirmed removed by the same grep and by
`tests/test_gui_smoke.py::test_baseline_removal_surface_is_absent` /
`test_baseline_module_is_removed`.
