# JOSS Readiness

An honest, evidence-based check against JOSS's actual review criteria
(https://joss.readthedocs.io/en/latest/review_criteria.html) as of this
writing. This document states what is in place and what is not — it does
not claim submission-ready status, since four concrete pending actions
(below) remain outside this repository's control.

## In place

| Criterion | Status | Evidence |
|---|---|---|
| Open source license | ✅ | GPL-3.0-or-later, OSI-approved (`LICENSE`) |
| Repository | ✅ | Public GitHub repository with version control history |
| Research application | ✅ | Storage-ring ion identification for nuclear/atomic physics; accompanying EPJ A manuscript (`RionID-EPJA/main.tex`, a companion work not distributed in this repository — see `docs/PUBLICATION_TRACEABILITY.md`) describes the method and an experimental application |
| Installation instructions | ✅ | `README.md` (PyPI and source), verified to actually work (`pip install -e ".[dev]"`, live-tested) |
| Example usage | ✅ | `examples/quickstart.py`, `docs/REPRODUCIBILITY.md` — both use synthetic data, no restricted inputs needed |
| Functionality documentation | ✅ | `README.md` parameter reference, `docs/SCIENTIFIC_METHOD.md`, `docs/LEGACY_BEHAVIOUR.md`, mkdocs API docs via mkdocstrings |
| Automated tests | ✅ | 19-test public pytest suite (`tests/`), covering the correction physics, extracted mass calculations, GUI state after feature removal, and the two speed-fixed hot paths' output identity |
| Continuous integration | ✅ | `.github/workflows/ci.yml` — test matrix (Python 3.9-3.12), lint, type-check, build, wheel-install smoke test |
| Contribution guidelines | ✅ | `CONTRIBUTING.md` |
| Community guidelines / reporting issues | ✅ | `CONTRIBUTING.md`, `SECURITY.md`, `.github/ISSUE_TEMPLATE/` |
| Citation metadata | ⚠️ partial | `CITATION.cff` exists but has a known, unresolved DOI discrepancy against the README's Zenodo badge — see below |
| Statement of need | ✅ | README "Scope and non-goals"; manuscript abstract/introduction |
| Substantial scholarly effort | ✅ | Multi-year development history (git tags v4.0.0-v8.0.0+), an accompanying peer-reviewable physics manuscript, and a contribution to a published experimental result (manuscript §5, citing the 2024 PRL two-photon-decay measurement) |

## Pending, not done in this repository

- **The CITATION.cff / README Zenodo DOI mismatch is not resolved.**
  `README.md`'s badge references `10.5281/zenodo.8169341`; `CITATION.cff`
  references `10.5281/zenodo.8169342`. This needs a decision from the
  repository owner about which is correct (or whether they refer to
  different things, e.g. concept vs. version DOI, and both should be
  stated explicitly) — see `docs/PUBLICATION_TRACEABILITY.md`. Not
  resolved by this packaging work, per explicit instruction.
- **No tagged 9.0.0 release exists yet.** `CITATION.cff`'s `date-released`
  still reads `2023-07-20` (from an earlier version) rather than a real
  9.0.0 release date, deliberately — see `CONTRIBUTING.md`'s release
  checklist for what happens at actual release time (tag, mint the
  version DOI via the Zenodo-GitHub integration, then update
  `CITATION.cff`/`.zenodo.json`).
- **No PyPI publish has occurred for 9.0.0.** `.github/workflows/publish.yml`
  exists and is `workflow_dispatch`-gated (never runs automatically) —
  publishing requires a deliberate, manual trigger plus a PyPI-side
  Trusted Publishing configuration that is also not done here.
- **JOSS submission itself has not occurred.** This document supports a
  future submission decision; it does not assert one has been made.

## Known limitations (documented, not fixed)

- **The `python3 -m rionid`/`rionid <datafile>` CLI entry point cannot
  currently complete a simulation run with any arguments.** Two
  independent, pre-existing bugs — discovered while writing
  `docs/REPRODUCIBILITY.md` and independently reproduced — mean `-psim`
  is effectively mandatory despite `argparse` not marking it required,
  and no `--circumference` flag exists at all, so every reference-
  frequency mode crashes. Neither is fixed in this packaging work (fixing
  requires changing `src/rionid/`, out of scope for a documentation-only
  plan); both are documented with the exact exception text and file
  provenance in `docs/OPEN_SCIENTIFIC_QUESTIONS.md` items 5-6. The GUI path is
  unaffected. `examples/quickstart.py` and `docs/REPRODUCIBILITY.md` §3
  both route around this by driving the underlying simulation engine
  directly rather than through the broken CLI path — this is why
  "Example usage" is still checked ✅ above: the shipped examples
  genuinely work, even though the raw CLI does not.

## Explicitly not claimed

This document does not claim: any adoption/usage statistics, any
performance benchmark beyond what `docs/PERFORMANCE_BASELINE.md` actually
measured (on one machine, on synthetic data, and stated as such), any DOI
beyond the two that already exist and conflict (see above), or that JOSS
review would find this repository acceptable — that determination belongs
to JOSS's editors and reviewers, not to this document.
