# Scientific Method

This is a reader-facing summary of the physics RionID implements and where
it lives in the code — written for a storage-ring scientist who wants to
use or verify the software, not as an implementation audit (that's
`docs/PUBLICATION_TRACEABILITY.md`). It cross-references the accompanying
EPJ A manuscript (`RionID-EPJA/main.tex`) by section and equation.

## 1. Reference-anchored forward model

Given a reference ion's mass-to-charge ratio μ_r and revolution frequency
f_r, the first-order relative frequency displacement of a candidate ion i
with mass-to-charge ratio μ_i is (manuscript Eq. `first_order`, §2.1):

```
f_i^(0) = f_r · [1 − α_p · (μ_i − μ_r) / μ_r]
```

where α_p is the momentum-compaction factor. This is implemented in
`ImportData.srrf` (`src/rionid/core.py`, `_calculate_srrf`) as a
dimensionless ratio array, one entry per candidate. Ionic mass-to-charge
ratios themselves come from `rionid.masses.ionic_moq_u`, which implements
the atomic-to-ionic mass correction (removed-electron rest mass, corrected
by the change in electron binding energy — manuscript §2.1, lines
198-200) using a verbatim-copied electron-binding-energy reference table
(cited in `masses.py`'s module docstring).

## 2. Empirical residual correction

Known anchor lines often reveal a smooth, setting-dependent departure from
the first-order model. RionID represents this as a residual — not an
absolute frequency — correction (manuscript Eq. `rev_correction`, §2.3):

```
Δf(x) = A·x² + B·x + C,     x = f^(0)
f^corr = f^(0) + Δf(f^(0))
```

Coefficient order is `(A, B, C)` = (quadratic, linear, constant), matching
`numpy.polyval([A, B, C], x)` — `A` has units Hz⁻¹, `B` is dimensionless,
`C` has units Hz. This is the **single, canonical implementation** of the
correction: `ImportData._calculate_srrf` (`src/rionid/core.py`), applied
in revolution-frequency space *before* harmonic multiplication. There is
no other place in the codebase that reimplements or duplicates this
arithmetic — `tests/test_correction.py` is a committed regression test
protecting it exactly.

The correction is empirical and setting-specific (manuscript §5): it
absorbs whatever smooth mismatch the first-order optics model, reference
choice, and beam conditions leave behind. It is not a measurement of any
particular physical quantity, and coefficients from one ring
configuration are not transferable to another without re-validation.

## 3. Harmonic projection

A pickup near harmonic h observes `F_{i,h} = h · f_i` (manuscript Eq.
`harmonic`, §2.2). RionID applies the correction once, in
revolution-frequency space, then multiplies by the requested harmonic(s) —
`ImportData._simulated_data` (`src/rionid/core.py`). If you need to import
a correction that was fitted directly in harmonic-frequency space at some
harmonic h₁, the transform to revolution-frequency space (and to any other
harmonic h₂) is given by manuscript Eq. `coefficient_transfer`, §2.4 — this
transform is a documentation/import convenience described in the
manuscript; the software's own internal representation always stays in
revolution-frequency space and never needs it at runtime.

## 4. What this software does not compute

Two things the manuscript discusses are explanatory/motivating, not
computed features of the software:

- **The harmonic-overlap criterion** (manuscript Eqs. `overlap_bounds`/
  `overlap_count`, §2.2) explains *why* candidate generation must consider
  a broad set of harmonics, not just the one nearest a pickup's nominal
  frequency — but RionID does not auto-derive which harmonics to display
  from this criterion. You supply the harmonic list explicitly (`-hrm`).
- **Uncertainty propagation** for the polynomial correction (manuscript
  Eq. `covariance`, §5) is described in the manuscript as a validation
  requirement, but is not currently computed or displayed anywhere in the
  software — there is no coefficient-covariance input path today. This is
  a known gap, not a hidden or partial implementation; see
  `docs/PUBLICATION_TRACEABILITY.md` for the full accounting.

## 5. Expert-guided assignment, not automatic classification

RionID computes and displays candidate positions; deciding which candidate
explains an observed line is left to you. The manuscript's own analyst
checklist (§3.2) — companion isotopic/isobaric pattern members, adjacent-
harmonic consistency, cross-detector agreement, charge-state/production
plausibility, competing-harmonic explanations — is not automated by this
software in any release. See `docs/AUTOMATIC_PID_REMOVAL_MAP.md` for the
record of automatic-matching functionality that existed in earlier
development and was deliberately removed before this release.
