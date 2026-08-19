# JOSS-packaging sub-project: design

Status: approved design, pre-implementation. Follows Wave 2a
(`docs/superpowers/specs/2026-08-19-wave2a-depid-masses-speed-design.md`,
merged to `master`). Assumes that spec and `docs/LEGACY_BEHAVIOUR.md`,
`docs/PUBLICATION_TRACEABILITY.md`, `docs/PERFORMANCE_BASELINE.md`,
`docs/AUTOMATIC_PID_REMOVAL_MAP.md`, `REFACTORING_PLAN.md` as read.

## Why this sub-project, and what it deliberately doesn't do

This is additive infrastructure and documentation work — no physics, no
scientific behaviour, no existing public API changes beyond the version
number. It does **not** touch `src/rionid/core.py`'s `_calculate_srrf`, any
other correction/calibration logic, or anything Wave 2a already removed.
Two items explicitly deferred out of Wave 2a's final review land here as
planned: adding CI, and setting up `ruff`.

Explicitly **not** in this sub-project (confirmed with you, or already
decided elsewhere):
- Wave 2b (the `correction.py`/`calibration.py`/`analysis.py`/etc. module
  split) — separate, independent sub-project, no ordering dependency on
  this one.
- The `exceptions.py` error-handling hierarchy — deferred to Wave 2b per
  `REFACTORING_PLAN.md`'s "explicitly out of scope" section.
- Any actual PyPI publish, GitHub release, or Zenodo archive action — this
  sub-project prepares the *infrastructure* for those (a non-triggering
  workflow, accurate placeholder files) but does not execute them.
- Resolving the CITATION.cff / README Zenodo DOI mismatch — stays flagged,
  per your explicit instruction during Phase 0.
- The EPJ A manuscript review/revision and its own journal-fit assessment —
  separate sub-project.
- Full type-annotation coverage of the codebase — `pyright` gets configured
  and made to pass in its default (lenient, "basic") mode; retrofitting type
  hints across an untyped codebase is out of scope.

**Version target: 9.0.0** everywhere version appears (`pyproject.toml`,
`CITATION.cff`, `src/rionid/version.py`), reflecting Wave 2a's breaking API
removal (automatic-PID and baseline-removal constructor/CLI parameters no
longer exist). `CHANGELOG.md` starts fresh at 9.0.0; earlier history is not
reconstructed from git log.

## Group A — Packaging & metadata

**`pyproject.toml`**: migrate to a proper PEP 621 `[project]` table as the
source of truth (name, version, description, readme, license, authors,
keywords, classifiers, `requires-python`, `dependencies`,
`[project.optional-dependencies]` with a `dev` group replacing
`[tool.poetry.group.dev.dependencies]`, and a `docs` group replacing
`[tool.poetry.group.docs]`). `[tool.poetry]` shrinks to just
`packages = [{include = "rionid", from = "src"}]` (still needed by
`poetry-core` for the src layout under `[build-system]`). `[project.urls]`
and `[project.scripts]` replace their `[tool.poetry.*]` equivalents.
**This is what actually makes `pip install -e ".[dev]"` work** — today it
can't, because no `[project.optional-dependencies]` exists at all under the
legacy pure-Poetry format.

Drop `black` from dev dependencies; the mega-brief's own Phase 5 checklist
names `ruff format --check .`, not black — one formatter, not two competing
opinions. Add `ruff`, `pyright`, `pre-commit`, `pip-audit`, `build`,
`twine` to the `dev` group (all needed for the CI/verification steps in
Group B and for the eventual Phase 5 check).

Acceptance for this specific change is testable, not just asserted: after
the edit, `python -m pip install -e ".[dev]"` in a scratch venv, `import
rionid`, the `rionid` console entry point resolving, and `python -m build`
producing an sdist+wheel that `twine check dist/*` accepts, must all
actually work in this sandbox (poetry-core is a PEP 517 backend; none of
this requires the `poetry` CLI itself, which isn't installed here).

**`CITATION.cff`**: version bumped to 9.0.0. DOI line untouched (mismatch
stays flagged, not resolved, per your instruction).

**`.zenodo.json`** (new): references the concept DOI as it appears in the
README today (`10.5281/zenodo.8169341`) via `related_identifiers` with
relation `isVersionOf` — doesn't introduce a third, different number, and
doesn't fabricate a version-specific DOI (Zenodo mints that automatically
on the next real release).

**`LICENSE`**: already GPL-3.0-or-later, OSI-approved. No action.

## Group B — Dev tooling & CI

**`ruff`** (`[tool.ruff]` in `pyproject.toml`): `target-version = "py39"`
(matches the declared support floor), line length 100, `select = ["E",
"F", "W", "I"]` (pycodestyle + pyflakes + isort) — ruff's own recommended
starting set, not maximum strictness. The handful of pre-existing
unused-import findings the final review already named (`core.py`'s
`append`/`peak_widths`, `gui/inputs.py`'s `QFont`/`time`,
`gui/plot.py`'s unused imports) get fixed as part of introducing the tool,
same spirit as Wave 2a's dead-code removal — not a broader stylistic
rewrite of untouched code.

**`pyright`** (`pyrightconfig.json` or `[tool.pyright]`):
`typeCheckingMode = "basic"` (pyright's lenient default — flags real
errors, not "missing type hints" as failures). If `pyright` surfaces
genuine errors (wrong arg counts, undefined names it can resolve, etc.),
fix them; if the only output is expected "no type information available"
noise in basic mode on untyped code, that's an acceptable passing state —
do not add type annotations project-wide to chase a stricter mode.

**`.pre-commit-config.yaml`**: `ruff` (lint + format) hooks, plus
`pre-commit-hooks`' `trailing-whitespace`, `end-of-file-fixer`,
`check-added-large-files`. Minimal, not exhaustive.

**`.github/workflows/ci.yml`** (new): matrix over Python 3.9–3.12 (the
declared `requires-python` range). Steps: checkout, setup-python, `pip
install -e ".[dev]"`, `ruff check .`, `ruff format --check .`, `pyright`,
then tests. Tests need `QT_QPA_PLATFORM=offscreen` (no display on CI
runners — the same headless approach validated throughout Wave 2a) and a
populated `~/.ame/` (all tests transitively call `get_ame_data()`, which
downloads from IAEA on a cold machine). Cache `~/.ame/` via
`actions/cache` keyed on a fixed key (the AME2020/mass table URL is
static) so only the first CI run per cache generation needs real network
access to IAEA; this is the fix for the exact gap the final review named.
A second job builds the sdist/wheel (`python -m build`, `twine check
dist/*`) and installs the built wheel into a fresh venv to smoke-test
`import rionid` and the `rionid` console entry point resolving (not a full
GUI launch — no display in CI — but enough to catch a broken package
build).

**Dependency/secret safety**: a `pip-audit` step in the CI workflow
(scans installed dependencies for known vulnerabilities — read-only,
doesn't touch the network beyond the audit database). `.github/dependabot.yml`
for automated dependency-update PRs (standard GitHub-native feature;
creates PRs, never pushes or merges anything itself).

**`.github/workflows/publish.yml`** (new, **non-triggering**): PyPI
Trusted Publishing (OIDC, no stored API token) via `pypa/gh-action-pypi-publish`,
gated on `workflow_dispatch` only — never runs on push or tag, requires a
human to manually invoke it from the GitHub Actions UI, and even then
requires the maintainer to have separately configured the trust
relationship on PyPI's own dashboard (an out-of-band action this sub-project
does not take and cannot take). Present as prepared, inert infrastructure.

## Group C — Documentation

**`README.md`** restructured with explicit sections: Purpose; **Scope &
Non-Goals** (states plainly that automatic/autonomous ion classification is
not part of this release, matching the manuscript's own claim); Installation;
Quick Start; Parameter Reference (drawing from `docs/LEGACY_BEHAVIOUR.md`'s
stage-by-stage parameter documentation); Supported Formats; Troubleshooting
(no-network/no-AME-cache case, the Qt-binding-conflict class of error
observed during Wave 2a generalised into a tip); Citation (points at
`CITATION.cff`); Limitations (drawn from the manuscript's own stated limits
and `docs/PERFORMANCE_BASELINE.md`'s open items); License; Acknowledgements
(kept as-is).

**`docs/SCIENTIFIC_METHOD.md`** (new): a reader-facing (not
audit-facing — distinct purpose from `PUBLICATION_TRACEABILITY.md`)
explanation of the physics — first-order frequency model, the polynomial
residual correction, the harmonic coefficient transform — with direct
cross-references to EPJ A section/equation numbers (via the existing
`PUBLICATION_TRACEABILITY.md` mapping) and to the code
(`core.py::_calculate_srrf`, `masses.py`'s ionic-mass correction). Written
against the *current* manuscript text, per your confirmation that the
physics/equations are stable even though the manuscript's prose will be
revised later.

**`docs/REPRODUCIBILITY.md`** (new): exact, runnable commands using the
now-existing public synthetic fixtures
(`tests/fixtures/synthetic_spectrum.py`) — expands
`PUBLICATION_TRACEABILITY.md`'s reproduction-commands section (written
during Phase 0, before the fixtures existed) with commands that actually
work today: running the golden correction test directly, generating a
synthetic spectrum and running the CLI against it, expected output.

**`examples/`** (new): one minimal, runnable example script using
`make_synthetic_spectrum` to generate data and `import_controller` (or the
CLI) to process it — clearly synthetic, no real experimental data,
matching the "minimal public example" requirement.

**`CHANGELOG.md`** (new, Keep a Changelog format): starts at 9.0.0 with a
detailed entry covering what Wave 2a removed (automatic PID, baseline
subtraction), changed (the `masses.py` extraction, the performance fixes),
and added (the public test suite). No attempt to backfill 4.0.0–8.0.0 from
git log.

**`CONTRIBUTING.md`** (new): dev setup (`pip install -e ".[dev]"`),
running tests/lint/type-check locally, PR process, and — folded in here
rather than as separate top-level files — a semantic-versioning policy
statement and a release checklist.

**`SECURITY.md`** (new): standard vulnerability-reporting template,
contact via the maintainer email already public in `CITATION.cff`/
`pyproject.toml`.

**`.github/ISSUE_TEMPLATE/bug_report.md`, `feature_request.md`,
`.github/PULL_REQUEST_TEMPLATE.md`** (new): standard, short templates.

## Group D — Synthesis

**`docs/JOSS_READINESS.md`** (new, written last, after A–C exist): a
checklist-style, evidence-based assessment against JOSS's actual review
criteria (license, research application, installation, examples, API/GUI
docs, tests, CI, contribution path, citation) plus an honest "not yet
done" list — tagging a release, minting the version DOI, resolving the
CITATION.cff DOI mismatch, and actually submitting are named as pending
maintainer actions, not claimed as complete. No invented adoption,
benchmark, or "submission-ready" status claims.

## Acceptance criteria

- [ ] `pip install -e ".[dev]"` works in a fresh venv in this sandbox.
- [ ] `python -m build && twine check dist/*` passes.
- [ ] `ruff check src/` and `ruff format --check src/` pass (scoped to
      `src/`, matching exactly what was cleaned -- see the plan's Task 5
      amendment for the full rationale).
- [ ] `pyright` runs and passes in basic mode (no genuine errors).
- [ ] `pytest` (the existing 19-test suite) still passes unmodified —
      this sub-project does not touch `tests/` or any file `tests/`
      exercises beyond `pyproject.toml`'s own metadata/tool config.
- [ ] CI workflow is present and its steps are individually verified
      runnable locally (GitHub Actions itself isn't invoked from this
      sandbox — no push happens — but every step it runs must have been
      shown to work standalone).
- [ ] The publish workflow exists, is `workflow_dispatch`-gated, and is
      never triggered by anything in this sub-project.
- [ ] No DOI, benchmark, adoption, or performance claim is fabricated
      anywhere in the new documentation.
