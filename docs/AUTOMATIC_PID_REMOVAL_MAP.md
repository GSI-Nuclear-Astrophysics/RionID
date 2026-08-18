# Automatic-PID Removal Map (Phase 0 Audit)

Status: **audit only — nothing has been removed yet**. This is the
inventory required before any deletion, per the task brief: every
PID-related file/function/UI item/config/test/dependency/claim, classified
as **Remove**, **Retain (rationale)**, or **Decision needed**. Removal
itself is a separate, later sub-project that will re-verify each item
against the passing regression suite before deleting it.

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

| # | Item | Location | Why |
|---|---|---|---|
| 1 | `setup_quick_pid()` | `gui/inputs.py:261-297` | Builds the "Quick PID" `QGroupBox` (alpha range, freq range, "Run Quick PID" button). |
| 2 | `quick_pid_script()` | `gui/inputs.py:422-537` | The scanning algorithm itself: nested loop over experimental peaks × α_p grid, calls `compute_matches`, silently overwrites `value_edit`/`alphap_edit`, emits the result as if user-confirmed. |
| 3 | `alphap_min_edit`, `alphap_max_edit`, `alphap_step_edit`, `fref_min_edit`, `fref_max_edit`, `threshold_edit` widgets, and their `load_parameters`/`save_parameters` keys (`alphap_min`, `alphap_max`, `alphap_step`, `fref_min`, `fref_max`, `threshold`) | `gui/inputs.py:88-93,119-124,231-236,265-290` | Config surface that exists only to feed #2. |
| 4 | `ImportData.compute_matches()` | `core.py:238-286` | Nearest-neighbour peak↔simulated-line χ² matcher; its only caller is #2. Also the function that overwrites `self.highlight_ions` with automatically matched species (`core.py:284`) — the clearest violation of "user provides desired ions; software calculates and displays them." |
| 5 | `ImportData.chi2`, `ImportData.match_count` attributes | `core.py:94-95` | State written only by #4. |
| 6 | `ImportData.save_matched_result()` | `core.py:288-295` | Exports the automatically-matched ion list. Confirmed **dead** even today (no caller anywhere in the GUI), but conceptually part of the automatic-matching feature — remove rather than wire up. |
| 7 | `matched_result_edit` field | `gui/inputs.py:87,118,254,256` | Loaded/saved to `parameters_cache.toml` but never connected to any action — a "half-supported flag" for #6. |
| 8 | `overlay_sim_signal` (`pyqtSignal`) + `MainWindow.overlay_simulation()` | `gui/inputs.py:30`, `gui/app.py:75,95-110` | Docstring: "used primarily during the 'Quick PID' scan to show visual feedback... without reloading the heavy experimental data every frame." Once #2 is gone, nothing emits this signal. The normal full-refresh path (`visualization_signal` → `MainWindow.update_visualization`) is unaffected and stays. |
| 9 | `RionID_GUI._stop_quick_pid`, `onPlotClicked()`, and the `MainWindow` wiring `visualization_widget.plotClicked.connect(rion_input.onPlotClicked)` | `gui/inputs.py:37,336-338,487,494,500`, `gui/app.py:78` | Stop-flag for #2's loop, set by clicking the plot. **Do not** remove `CreatePyGUI.plotClicked` itself (`gui/plot.py:37,58`) — it also drives the legitimate manual "Pick" coordinate feature (see Retain list). Only this specific subscriber goes. |
| 10 | README.md / `docs/index.md` "Automated Matching: Includes Quick PID logic to scan α_p and Reference Frequency to find the best match (χ² minimization)" bullet | `README.md`, `docs/index.md` | Unsupported-capability claim for this release; contradicts the manuscript's own "not an autonomous classifier" statement. |

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

| # | Item | Location | Resolution |
|---|---|---|---|
| 1 | `Particle.identify_range`, `identify_search_range`, `identify_region`, `identify_search_region`, `get_unknown_moq_from_freq`, `get_unknown_rev_freq_from_moq`, `get_nuclides_freqs` | `external/barion/particle.py:392-474` | **Remove.** You are also an owner of the upstream `barion` library, so vendoring/re-vendoring risk does not apply the way it would for an arbitrary third party. Decision: **drop the `external/barion` vendored subpackage as a dependency entirely**; extract only the routines RionID actually calls (ionic-mass/m-q calculation including the electron-binding correction, AME/NUBASE mass-table loading and unit constants, and a minimal ring-circumference holder) into RionID's own source tree. `identify_*`/`get_unknown_*`/`get_nuclides_freqs` are not among the extracted routines and are dropped. |
| 2 | `Particle.get_all_in_all`, `get_isotopes`, `get_isotones`, `get_isobars`, `Ring.get_ring_dict()` | `external/barion/particle.py:83-171`, `external/barion/ring.py:29-57` | **Remove**, same resolution as #1 — not part of the extracted subset RionID uses. |
| 3 | `handle_tiqnpz_data`, `handle_prerionidnpz_data` | `io.py:72-101,124-141` | **Remove.** Confirmed dead code; no re-wiring requested. |
| 4 | `CollapsibleGroupBox` | `gui/inputs.py:546-574` | **Remove**, unused, unrelated to PID. |

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
removal is needed or possible from this change alone.**

## Tests

No automatic-PID tests exist (no tests exist at all yet — see
`PUBLICATION_TRACEABILITY.md`). Nothing to remove here; the upcoming
regression suite (Phase 1) should simply never test `compute_matches` /
`quick_pid_script`, since they will not exist.

## Documentation/claims to update once removal lands

- `README.md` "Features" bullet (item 10 above).
- `docs/index.md` (identical bullet, generated from the same README content).
- `RionID-EPJA/main.tex` already correctly claims *no* automatic
  classification — no manuscript text changes are needed for this reason
  specifically (Phase 4 will still review the manuscript for other reasons).
- `parameters_cache.toml` (local, untracked) will simply stop being written
  the removed keys once the code is gone; no action needed.
