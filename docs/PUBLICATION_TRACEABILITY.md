# Publication Traceability (Phase 0 Audit)

Status: **audit only — no source or manuscript changes**. Maps every
equation, figure, table, and load-bearing claim in
`RionID-EPJA/main.tex` (766 lines, compiled clean to a 7-page PDF with
`latexmk` during this audit — see "TeX/PDF reconciliation" below) to the
code that implements it, its current test coverage, and a reproduction
command. Precedence follows the order given in the task brief: explicit
instructions → manuscript source → compiled PDF → verified code → historical
notebooks.

## TeX/PDF reconciliation

`RionID-EPJA/main.tex` was compiled locally with `latexmk -pdf` (LaTeX
Modern fonts, `svjour3`-adjacent `article` class, `10pt,twocolumn`). Result:
7 pages, no undefined citations (all 10 `\cite` keys resolve against the 11
manually-typed `\bibitem` entries in `main.tex:692-764`), no missing
references. Two reconciliation notes:

1. `RionID-EPJA/bibliography.bib` (17 KB, present alongside `main.tex`) is
   **not referenced by the document** — the manuscript uses a manual
   `thebibliography` environment (`main.tex:692`), not `\bibliography{}`.
   The `.bib` file appears to be a separate/orphaned reference-manager
   export. Not a defect, but worth resolving before submission (keep as a
   Zotero/BibTeX source of truth and regenerate `thebibliography` from it,
   or delete it) — **flagged for your decision, not changed**.
2. The manuscript's own header comment (`main.tex:1-14`) already lists five
   "REQUIRED AUTHOR ACTIONS BEFORE SUBMISSION" (version DOI, example-data
   DOI, author/CRediT/funding confirmation, numerical calibration
   diagnostics from real data, move into the Springer/EPJ A template) and
   defines `\todoentry{}`-flagged placeholders `\paperreleasedoi` and
   `\exampledatadoi`. These are the author's own tracked open items, not
   discoveries of this audit — repeated here only so Phase 4 work treats
   them as authoritative.

## Equation-by-equation map

| Manuscript | Code | Status |
|---|---|---|
| Eq. (1) `first_order_full` (general first-order model, incl. velocity term) | Not implemented — code uses the reduced form only | **Expected**: text explicitly presents Eq. (1) as the general case and Eq. (2) as what RionID actually uses ("under sufficiently cooled conditions... RionID uses Eq. (2)"). Consistent. |
| Eq. (2) `first_order` (`f_i^(0)`) | `ImportData.srrf`, `core.py:352-353` | **Verified** — direct transcription, confirmed by reading. |
| Eq. (3) `harmonic` (`F=hf`) | `ImportData._simulated_data`, `core.py:393` | **Verified**. |
| Eqs. (4)-(7) `rev_interval`/`harmonic_interval`/`overlap_set`/`overlap_bounds` (geometric harmonic-overlap criterion) | **Not implemented anywhere in the software.** | The manuscript presents this as an *explanatory/motivating* criterion (why candidate generation must consider multiple harmonics) illustrated in Fig. 2 (`count_overlaps.png`), not as a computed, exposed feature. `-hrm/--harmonics` (CLI) / "Harmonics:" (GUI) is a plain user-supplied list — RionID does not auto-derive it from Eqs. (4)-(7). **This is a scope clarification for `SCIENTIFIC_METHOD.md`, not a bug**: the manuscript should say plainly that the overlap criterion motivates *why* users must supply a broad-enough harmonic list, not that the software computes it for them. |
| Eq. (8) `overlap_count` (`N_h`) | Same as above — not implemented. | Numeric example (`N_5=0, N_125=8, N_210=15`, `main.tex:275-276`) and Fig. 2 appear to have been produced by an external/one-off calculation, not the shipped package. **Reproducibility gap**: flag for Phase 3/4 — either add a small `analysis.py` helper that computes Eq. (8) so Fig. 2 is regenerable from the package, or state explicitly in the paper that Fig. 2 was produced by auxiliary code outside the release. **Needs your decision.** |
| Eq. (9) `rev_correction` (`Δf(x)=Ax²+Bx+C`) | `ImportData._calculate_srrf`, `core.py:354-356` | **Verified numerically** in this audit (see below). |
| Eq. (10) `corrected_rev` (`f^corr = f^(0)+Δf`) | Same | **Verified**. |
| Eq. (11)-(13) coefficient-order/units convention | `core.py:355` (`polyval(array(correct), ...)`, `-c a0 a1 a2` CLI help, GUI placeholder "a0\*x\*\*2 + a1\*x + a2") | **Verified** — `numpy.polyval` convention matches the paper's stated `(A,B,C)=(quadratic,linear,constant)` exactly. |
| Eqs. (14)-(19) `harmonic_scaling`/`coefficient_h`/`coefficient_transfer`/`scaling_invariant`/`h214_numbers` | **Not implemented as a callable function** — the code always works in revolution-frequency space (harmonic=1) internally and never needs the harmonic→harmonic coefficient transform at runtime. | This is a **documentation/worked-example** contribution of the paper (for users importing legacy per-harmonic coefficients), not a code path. **Verified numerically anyway** (see below) because it is the strongest available check that the code's correction convention is self-consistent with the paper's stated math. |
| Eq. (20) `h127_numbers`/(21) `h127_corrected` | — | Reproduced exactly by direct arithmetic in this audit; see next section. |
| Eq. (22) `covariance` (`σ²_Δf(x) = v(x)ᵀ Σ_β v(x)`) | **Not implemented anywhere.** | **Significant gap.** No function in `core.py` (or elsewhere) computes or exposes correction-coefficient covariance or a propagated uncertainty. The product brief explicitly requires the app to "show essential fit/correction information, units, **uncertainties**, warnings, and provenance" — today's GUI shows none of that for the polynomial correction. This needs a design decision in Phase 2 (what covariance input format, what display) — **flagged, not silently added**. |

## Numerical verification performed in this audit

Using the manuscript's own harmonic-214 numbers (`main.tex:372-377`) and the
transform of Eqs. (14) and (17):

```
a214, b214, c214 = 1.19262945e-9, -9.85167644e-1, 2.03447706e8   # paper
A = 214*a214 = 2.55222702e-07   B = b214 = -9.85167644e-01   C = c214/214 = 9.50690215e+05
  -> matches main.tex:380-382 to all printed digits.
a127 = A/127 = 2.00962758e-09   b127 = B   c127 = C*127 = 1.20737657e+08
  -> matches main.tex:386-389 to all printed digits.
```
Applying the code's actual operation order — correct in revolution-frequency
space, then multiply by harmonic (`core.py:354-356` then `core.py:393`) —
against the direct harmonic-127 polynomial of Eq. (21) gives a residual of
`0.000e+00 Hz` at an arbitrary test frequency. **Conclusion: the shipped
`_calculate_srrf` implementation is numerically consistent with the
manuscript's correction formalism, including the harmonic-transform
identity.** Recommend this exact check become a committed regression test
(`tests/test_correction.py`) in Phase 1 — it is the single highest-value
golden-number test available, and it currently exists only as an
audit script, not as project-owned, CI-run code.

## Section-level map

| Manuscript section | Code / UI | Status |
|---|---|---|
| §3.1 "five groups of information" (spectrum, candidate list, reference ion+α_p, harmonics, correction) | CLI: `-r/-ap/-hrm` + one of `-b/-ke/-gam/-f` + `-psim` + `-c`; GUI: matching fields in `gui/inputs.py` | **Verified**, 1:1. |
| §3.1 deterministic 6-step sequence (`main.tex:430-440`) | `_calculate_moqs` → `_calculate_srrf` → `_simulated_data` → (display filter) → render | **Verified**, order matches. |
| §3.2 "Expert-guided assignment" (5-point analyst checklist) | **Intentionally not automated** — GUI supplies overlay/cross-detector viewing only | **Consistent** with manuscript framing; see next entry for the one place this consistency currently breaks. |
| Abstract / §5 Conclusions: *"not... an autonomous classifier"*, *"expert-guided assignment... [not] a replacement for a dedicated mass calibration"* | **Contradicted by shipped code**: `gui/inputs.py.quick_pid_script` (`core.py.compute_matches`) implements exactly an automatic χ²-minimising optimiser over (α_p, reference frequency) that silently overwrites user-entered values and the `highlight_ions` selection. README/`docs/index.md` explicitly advertise it ("Automated Matching... χ² minimization"). | **This is the central finding motivating the automatic-PID removal task** — see `docs/AUTOMATIC_PID_REMOVAL_MAP.md`. After removal, the claim becomes true of the shipped software; today it is not. |
| §4 "Reproducible use case" / Program Summary `\exampledatadoi` | **Does not exist yet** — no `tests/`, `examples/`, or archived fixture set in the repository | Major gap for both JOSS and EPJ A; central task for the regression-fixture sub-project (Phase 1/3). |
| Program Summary "Version documented: 8.0.0" | `pyproject.toml: version = "8.0.0"`, `src/rionid/version.py` | **Verified**, matches exactly. |
| Program Summary / Code availability: repository URL, licence GPL-3.0-or-later | `pyproject.toml`, `LICENSE` | **Verified**, matches (GitHub org/repo casing differs cosmetically: `GSI-Nuclear-Astrophysics/RionID` in the paper vs `.../rionid` in `pyproject.toml`; GitHub URLs are case-insensitive, not a functional issue). |
| Persistent Zenodo concept DOI `10.5281/zenodo.8169341` (`main.tex:114`, README badge) | `CITATION.cff: doi: 10.5281/zenodo.8169342` | **MISMATCH** — differs in the final digit. One of the two is wrong, or they intentionally refer to different Zenodo records (concept vs. version DOI) and the CFF should say so explicitly. **Flagged for your decision; not silently resolved.** |
| §5 validation requirements (algebraic harmonic-invariant tests, regression tests vs. archived line list, out-of-fit held-out anchor, CI across supported Python versions) | **None of this exists yet** | Directly defines the Phase 1 regression-test scope: (1) the harmonic-invariant algebraic test (now verified manually, needs to become code), (2) a fixed-candidate-input regression test against an archived expected line list (needs a public fixture), (3) CI matrix across the declared Python 3.9-3.12 range (current CI only deploys docs, `.github/workflows/deploy.yml`). |
| §4 "same workflow contributed to... two-photon-decay experiment" citation to `FreirePRL2024` | External, not verifiable from this repository | Out of scope for code audit; noted for completeness. |

## Figures and tables

| Figure | Source data | Reproducible from repo today? |
|---|---|---|
| Fig. 1 `1d-full-labels.png` (broadband E143 spectrum, h=123-129 triplet) | Real E143 data, not in repo | No — requires the archived reduced dataset promised at `\exampledatadoi`. |
| Fig. 2 `count_overlaps.png` (Eq. 8 staircase) | Not computed by the package (see Eq. 8 row above) | No — and not reproducible from the package even with data, since Eq. (8) isn't implemented. |
| Fig. 3 (a)/(b) `245vs410.png`/`410vs245.png` (dual-pickup contaminant) | Real E143 data | No, same reason as Fig. 1. |
| Fig. 4 `correct_identification72Ge.png` (calibration residuals) | Real E143 anchor data; caption itself says "must be reported from the source data before submission" | No — and the caption already flags this as an open item, consistent with this audit's findings. |

No tables (numeric results appear inline in Eqs. 15-21, already verified
above).

## Test coverage today

**None.** There is no `tests/` directory in the repository at all. Every
"Verified" status above reflects manual reading and the ad hoc numeric
checks performed during this audit, not committed, CI-run tests. Building
the regression suite that encodes these checks (starting with the
correction-formula and harmonic-invariant tests, which are already fully
specified by the manuscript) is the first concrete deliverable of the next
sub-project.

## Reproduction commands (current, pre-refactor CLI)

```bash
# Deterministic forward simulation + correction, no GUI, writes simulation_result.out
python -m rionid <datafile>.npz -r 72Ge+35 -ap 0.189 -f 1930000 -hrm 127 -c 2.55e-7 -0.985 9.5e5

# Same, with the GUI shown
python -m rionid <datafile>.npz -r 72Ge+35 -ap 0.189 -f 1930000 -hrm 123 124 125 126 127 128 129 -s
```
These are the CLI's actual current argument names (`-r/-ap/-f/-hrm/-c/-s`),
confirmed against `__main__.py:24-46`; no `<datafile>.npz` fixture ships
with the repository today.
