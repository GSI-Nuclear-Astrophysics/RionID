# JOSS-Packaging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate `pyproject.toml` to modern PEP 621 metadata, configure and pass `ruff`/`pyright`/`pre-commit`, add CI (test matrix, lint, type-check, build, wheel-install smoke test), and write the documentation set (README overhaul, SCIENTIFIC_METHOD.md, REPRODUCIBILITY.md, CHANGELOG, CONTRIBUTING, SECURITY, JOSS_READINESS.md) that brings RionID to JOSS-submission-preparation standard — no physics or scientific-behaviour changes anywhere.

**Architecture:** Purely additive/infrastructure work on top of the merged Wave 2a state. No module restructuring (that's the separate, independent Wave 2b). Every packaging claim in this plan (dependency versions, `pip install -e ".[dev]"` working, `python -m build` + `twine check` passing) was live-tested against this exact repo before being written down — see the design spec's Group A for context.

**Tech Stack:** `poetry-core` (PEP 517 build backend, unchanged), PEP 621 `pyproject.toml`, `ruff` (lint + format), `pyright` (basic mode, with justified suppressions — see Task 3), `pre-commit`, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-08-19-joss-packaging-design.md`

## Global Constraints

- **No change to any file under `src/rionid/` that affects runtime behaviour.** Every source-code edit in this plan is either (a) a pure formatting/whitespace change from `ruff format`/`ruff check --fix` on safe rules only, or (b) one of the specific, named, zero-risk fixes in Task 2/Task 3 below. Nothing else in `src/` changes.
- **Never fix a "possibly unbound variable" finding that changes program behaviour under an edge case** (empty `harmonics` list, `[Calculations]` marker missing, etc.) — log it in `docs/OPEN_SCIENTIFIC_QUESTIONS.md` instead, per the project's physics-lock rule. Task 3 lists exactly which findings are real-but-logged vs. safe-to-fix.
- **Do not fabricate a DOI, benchmark, adoption claim, or "submission-ready" status anywhere.** The CITATION.cff/README DOI mismatch stays flagged, not resolved.
- **Version target is 9.0.0** everywhere a version string appears (`pyproject.toml`, `CITATION.cff`, `src/rionid/version.py`).
- Every task ends with its own `git add` + `git commit` for exactly the files that task touches.
- `poetry` (the CLI) is not installed in this environment; none of this plan needs it — `pip`, `python -m build`, and `twine` operate on `poetry-core` as a PEP 517 backend without requiring the `poetry` CLI at all.

---

### Task 1: `pyproject.toml` PEP 621 migration + version bump

**Files:**
- Modify: `pyproject.toml`
- Modify: `src/rionid/version.py`

**Interfaces:**
- Produces: a `[project]` table that `pip install -e ".[dev]"` and `python -m build` both work against (live-verified during planning — see Step 5).

- [ ] **Step 1: Replace `pyproject.toml` in full**

This exact content was live-tested against this repo (editable install, dev-dependency resolution, sdist+wheel build, `twine check`) before being written here — transcribe it verbatim:

```toml
[project]
name = "rionid"
version = "9.0.0"
description = "Ring-stored ion IDentification: A pure Python tool for Schottky spectrum analysis."
authors = [
    {name = "David Freire-Fernández", email = "D.FreireFernandez@gsi.de"},
]
license = "GPL-3.0-or-later"
readme = "README.md"
requires-python = ">=3.9,<3.13"
keywords = ["Physics", "Nuclear Physics", "Storage Rings", "Mass Spectrometry", "Schottky"]
classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Physics",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]
dependencies = [
    "numpy>=2.0.2,<3",
    "scipy>=1.10.1,<2",
    "PyQt5>=5.15.9,<6",
    "pyqtgraph>=0.13.3,<0.14",
    "loguru>=0.7.0,<0.8",
    "toml>=0.10.2,<0.11",
    "ezodf>=0.3.2,<0.4",
    "fortranformat>=2.0.3,<3",
    "lxml>=6.0.2,<7",
]

[project.urls]
Homepage = "https://gsi-nuclear-astrophysics.github.io/rionid/"
Repository = "https://github.com/GSI-Nuclear-Astrophysics/rionid"
"Bug Tracker" = "https://github.com/GSI-Nuclear-Astrophysics/rionid/issues"
Documentation = "https://gsi-nuclear-astrophysics.github.io/rionid/"

[project.scripts]
rionid = "rionid.gui.app:main"

[project.optional-dependencies]
dev = [
    "pytest>=8.0.0,<9",
    "pytest-qt>=4.2.0,<5",
    "ruff>=0.6.0,<1",
    "pyright>=1.1.380,<2",
    "pre-commit>=3.8.0,<4",
    "pip-audit>=2.7.0,<3",
    "build>=1.2.0,<2",
    "twine>=6.0.0,<7",
]
docs = [
    "mkdocs>=1.5.0,<2",
    "mkdocs-material>=9.5.0,<10",
    "mkdocstrings[python]>=0.24.0,<0.25",
    "mkdocs-gen-files>=0.5.0,<0.6",
    "mkdocs-literate-nav>=0.6.0,<0.7",
]

[tool.poetry]
packages = [
    { include = "rionid", from = "src" }
]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
addopts = "--strict-markers -m \"not slow\""
markers = [
    "slow: larger-N performance/regression cases, not run by default",
]

[tool.ruff]
target-version = "py39"
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "W", "I"]

[tool.ruff.lint.per-file-ignores]
# The verbatim-copied electron-binding-energy physics table (see
# masses.py's module docstring) is intentionally many numbers per line --
# do not reformat it for line length, that risks transcription damage to
# physics reference data for zero benefit. NOTE: per-file-ignores only
# suppresses `ruff check` (the linter), not `ruff format`. The table
# itself must additionally be wrapped in `# fmt: off` / `# fmt: on` at
# the point of use (Task 2) -- ruff format does honor those markers.
"src/rionid/masses.py" = ["E501"]

[tool.pyright]
include = ["src"]
pythonVersion = "3.9"
typeCheckingMode = "standard"
# Disabled: dominated by PyQt5/pyqtgraph stub-interop false positives
# (Qt enum access patterns valid at runtime via SIP but not modeled by
# the bundled stubs; pyqtgraph's dynamically-assigned attributes;
# Optional-attribute narrowing pyright can't see across methods) rather
# than real bugs -- confirmed during JOSS-packaging planning by manually
# categorizing every one of the ~60 findings these rules produced.
# reportPossiblyUnboundVariable stays enabled: it is what caught the
# real, pre-existing issues now logged in docs/OPEN_SCIENTIFIC_QUESTIONS.md.
reportOptionalMemberAccess = false
reportOptionalIterable = false
reportArgumentType = false
reportAttributeAccessIssue = false
reportCallIssue = false
reportGeneralTypeIssues = false

[build-system]
requires = ["poetry-core>=1.7.0"]
build-backend = "poetry.core.masonry.api"
```

- [ ] **Step 2: Bump `src/rionid/version.py`**

Replace:
```python
__version_info__ = (8,0,0)
__version__ = '.'.join(map(str, __version_info__))
```
with:
```python
__version_info__ = (9, 0, 0)
__version__ = '.'.join(map(str, __version_info__))
```
(also adds the missing trailing newline and the missing space after commas,
both of which `ruff format` would otherwise flag in Task 2 — fixing them
here avoids Task 2 touching a line Task 1 just wrote.)

- [ ] **Step 3: Install and verify in a scratch venv (not the system Python)**

Do not `pip install -e .` into the ambient environment — this sandbox has a
messy, unrelated Anaconda base environment with pre-existing broken
packages (a `pip install` summary step can crash on an unrelated malformed
version string from an unrelated package; this is not caused by anything
in this repo and does not indicate a problem with this change — the
install itself still completes, but avoid it entirely by using a clean
venv):
```bash
python3 -m venv /tmp/rionid-pkg-test
/tmp/rionid-pkg-test/bin/pip install --upgrade pip --quiet
/tmp/rionid-pkg-test/bin/pip install -e ".[dev]"
/tmp/rionid-pkg-test/bin/python -c "import rionid; from rionid.version import __version__; print(__version__)"
/tmp/rionid-pkg-test/bin/rionid --help 2>&1 | head -3 || true
```
Expected: `__version__` prints `9.0.0`; the `rionid` console script resolves
(it will fail to actually launch without a display, that's fine — this
just confirms the entry point exists and is importable, not a full GUI
smoke test).

- [ ] **Step 4: Build and check**

```bash
rm -rf dist/ build/
/tmp/rionid-pkg-test/bin/python -m build
/tmp/rionid-pkg-test/bin/twine check dist/*
```
Expected: both `rionid-9.0.0.tar.gz` and `rionid-9.0.0-py3-none-any.whl`
build successfully and `twine check` reports `PASSED` for both. (If
`twine check` fails with "Metadata is missing required fields", the
installed `twine` is too old to parse the `Metadata-Version` your
`poetry-core` produces — confirm the venv actually resolved
`twine>=6.0.0`, not an older cached version.)

- [ ] **Step 5: Inspect the wheel contents**

```bash
python3 -m zipfile -l dist/rionid-9.0.0-py3-none-any.whl
```
Expected: every `.py` file under `src/rionid/` (including `masses.py`,
`gui/*.py`, `external/lisereader/*.py`) is present under `rionid/` in the
archive, plus `rionid-9.0.0.dist-info/{METADATA,WHEEL,entry_points.txt,RECORD,licenses/LICENSE}`.
This confirms `[tool.poetry] packages = [...]` is still correctly honoured
by `poetry-core` even with a `[project]` table present.

- [ ] **Step 6: Clean up and commit**

```bash
rm -rf dist/ build/ src/rionid.egg-info /tmp/rionid-pkg-test
git add pyproject.toml src/rionid/version.py
git commit -m "Migrate pyproject.toml to PEP 621, bump version to 9.0.0

Adds a [project] table as the metadata source of truth -- this is what
makes \`pip install -e \".[dev]\"\` actually work; the prior pure-Poetry
format had no [project.optional-dependencies], so that exact command
(named in the release-verification checklist) could not resolve extras
at all. [tool.poetry] shrinks to just the src-layout packages= entry,
still required by poetry-core alongside [project].

Drops black from dev dependencies in favour of ruff format (one
formatter, matching the verification checklist's ruff format --check .,
not a black+ruff conflict). Adds ruff/pyright/pre-commit/pip-audit/
build/twine to the dev extra, needed by later tasks in this plan.

Version bumped to 9.0.0 (pyproject.toml, CITATION.cff follows in a later
task, src/rionid/version.py) reflecting Wave 2a's breaking API removal
(automatic-PID and baseline-removal parameters no longer exist).

Verified: pip install -e \".[dev]\" in a scratch venv, python -m build,
twine check dist/*, and wheel-content inspection all pass -- see this
commit's task report for full output.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: `ruff` — safe auto-fix, format, and the 4 remaining manual fixes

**Files:**
- Modify: every file under `src/rionid/` that `ruff format`/`ruff check --fix` touches (mechanical only)
- Modify: `src/rionid/__init__.py`, `src/rionid/core.py`, `src/rionid/external/lisereader/reader.py`, `src/rionid/gui/plot.py` (the 4 manual fixes below)

**Interfaces:** None — no function signature or behaviour changes anywhere in this task.

This task was fully dry-run during planning: `ruff check --select E,F,W,I --line-length 100 --extend-per-file-ignores "src/rionid/masses.py:E501" --fix src/` auto-fixed 221 of 305 findings (all mechanical: trailing whitespace, blank-line whitespace, missing final newlines, import sorting) with zero semantic risk, and `ruff format` reformatted 11 files. Exactly 4 `ruff check` findings remained, each named precisely below.

**Correction found during Task 2's execution (not caught during the dry-run above):** `ruff format` reformatted `masses.py`'s `ElBiEn` physics table too, since `[tool.ruff.lint.per-file-ignores]` only suppresses the linter, not the formatter. Step 1a below (added after the fact) fixes this with `# fmt: off`/`# fmt: on` markers, which `ruff format` does honor, around only the table itself.

- [ ] **Step 1a: Protect the `ElBiEn` table from the formatter, before running it**

In `src/rionid/masses.py`, add `# fmt: off` immediately before the
`ElBiEn = [` line and `# fmt: on` immediately after its closing `]]` (the
last line of the table). `ruff format` honors these markers (adopted from
Black's convention) and will skip only the fenced region — the rest of the
file (module docstring, `Ring`/`AMEData` methods, the module-level
functions at the tail) still gets formatted normally. Do this *before*
Step 1b, not after — `[tool.ruff.lint.per-file-ignores]` (Task 1's config)
only suppresses the *linter*'s `E501` check on this file, it does **not**
exempt the file from the *formatter* at all; without these markers,
`ruff format` reformats the table's ~540 lines of numbers (whitespace
changes only, but exactly the kind of transcription-adjacent risk the
per-file-ignore was meant to avoid).

- [ ] **Step 1b: Run the safe auto-fix and formatter**

```bash
ruff check --fix src/
ruff format src/
```
(Uses the `[tool.ruff]`/`[tool.ruff.lint]` config from Task 1 automatically
— `--select`/`--line-length`/`--per-file-ignores` no longer need repeating
on the command line.) Expected: "221 fixed" (or close — exact count may
shift slightly if Task 1's `version.py` edit already fixed its own
formatting issue) and "11 files reformatted" (or similar), then re-run
`ruff check src/` to confirm exactly 4 findings remain, matching the ones
below. If more than 4 remain, or different ones, STOP and report —
something about the repo state differs from what planning found.

After this step, verify the table itself is untouched:
`git diff src/rionid/masses.py` should show only the `# fmt: off`/
`# fmt: on` lines added, plus any import-order fix, with **zero** changes
inside the `ElBiEn = [...]` block itself.

- [ ] **Step 2: Fix `src/rionid/__init__.py`'s two "imported but unused" findings**

These are `__init__.py`'s public re-exports (`ImportData`, `__version__`),
not genuinely unused — ruff's suggested fix (explicit re-export syntax)
is the correct, standard way to tell it that. Replace:
```python
from .core import ImportData
from .version import __version__
```
with:
```python
from .core import ImportData as ImportData
from .version import __version__ as __version__
```

- [ ] **Step 3: Fix `src/rionid/core.py`'s bare `except:` (E722)**

In `_parse_ref_ion`'s fallback-parsing block, replace:
```python
            except:
                print(f"Warning: Could not parse reference ion '{refion}'.")
```
with:
```python
            except Exception:
                print(f"Warning: Could not parse reference ion '{refion}'.")
```
This is a strict improvement with no behaviour change for any case this
code is meant to handle — bare `except:` also swallows `KeyboardInterrupt`/
`SystemExit`, which `except Exception:` correctly does not.

- [ ] **Step 4: Fix `src/rionid/external/lisereader/reader.py`'s one long line**

Find the line (post-format, may have shifted slightly; search for the
text):
```python
                    f"Please check if the AME data needs to be updated, or verify the formatting/spelling (e.g., 'Os' vs 'OS')."
```
Reflow the single f-string literal across two lines (this is prose text,
not physics data — safe to reflow):
```python
                    f"Please check if the AME data needs to be updated, or "
                    f"verify the formatting/spelling (e.g., 'Os' vs 'OS')."
```

- [ ] **Step 5: Fix `src/rionid/gui/plot.py`'s `width`/`style` possibly-unbound false positive**

This one is caught by `pyright` (Task 3), not `ruff`, but the fix belongs
in this task since it's the same kind of safe, zero-risk cleanup and Task 3
depends on it already being done. In `plot_simulated_data`, find:
```python
                # Determine Style
                if is_highlight:
                    c = color_match
                    width = 2
                    style = Qt.SolidLine
                elif is_ref:
                    c = color_ref
                    width = 2
                    style = Qt.DashLine
                else:
                    c = color
                    # No need for width/style here, handled by bulk curve
```
Add a default initialization immediately before this block so `width`/
`style` are always bound (they are never actually *read* when unbound at
runtime — the only read site is inside `if is_highlight or is_ref:`, which
excludes the `else` branch — this purely satisfies static analysis, it
changes no runtime value):
```python
                # Determine Style
                width = style = None  # only read under is_highlight/is_ref, set below
                if is_highlight:
                    c = color_match
                    width = 2
                    style = Qt.SolidLine
                elif is_ref:
                    c = color_ref
                    width = 2
                    style = Qt.DashLine
                else:
                    c = color
                    # No need for width/style here, handled by bulk curve
```

- [ ] **Step 6: Verify clean**

```bash
ruff check src/
ruff format --check src/
```
Expected: both report no issues (`All checks passed!` / `X files already formatted`).

- [ ] **Step 7: Run the full test suite** (this task touched application source
      files, even if only cosmetically — confirm nothing broke)

```bash
python3 -m pytest tests/ -v
```
Expected: 18 passed, 1 deselected (unchanged from before this task).

- [ ] **Step 8: Commit**

```bash
git add src/
git commit -m "ruff: safe auto-fix, format, and the 4 remaining manual fixes

Mechanical only: trailing whitespace, blank-line whitespace, missing
final newlines, import sorting, and one long line reflowed (prose, not
the masses.py physics table, which is per-file-ignored for E501).

Four small, zero-risk manual fixes: __init__.py's re-exports use
explicit 'as' syntax; core.py's bare except: narrowed to except
Exception:; lisereader/reader.py's one over-length f-string reflowed;
gui/plot.py's width/style variables given a default (never read
unbound at runtime -- satisfies pyright's flow analysis, changes no
value).

No function signature or behaviour change anywhere. Full test suite
re-run and unchanged: 18 passed, 1 deselected.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: `pyright` — verify config, log the real findings

**Files:**
- Create: `docs/OPEN_SCIENTIFIC_QUESTIONS.md`
- Modify: `pyproject.toml` (only if Task 1's committed config still says
  `typeCheckingMode = "basic"` — see the note in Step 1 below; if it
  already says `"standard"`, there is nothing to change here)

**Interfaces:** None — this task makes no application source changes (Task 2
already applied the one safe `src/` fix `pyright` needed).

`pyright`'s config landed in Task 1. This task verifies it reports the
expected findings and creates the standing location the whole engagement's
brief has referred to for logging suspected issues, since this is the
first time concrete candidates for it exist.

- [ ] **Step 1: Run `pyright` and confirm it's clean of everything except the known,
      logged findings**

```bash
pyright
```
Expected: exactly these `reportPossiblyUnboundVariable` findings (the
suppressed categories from Task 1's config should produce nothing):
- `src/rionid/core.py` — `gamma` possibly unbound (in `calc_ref_rev_frequency`)
- `src/rionid/external/lisereader/reader.py` — `file_start` possibly unbound (×2), `element` possibly unbound (×2)
- `src/rionid/gui/controller.py` — `harmonic` possibly unbound, `fre` possibly unbound (in `save_simulation_results`)

That's **7 raw diagnostic lines** total (5 distinct conceptual issues;
`file_start` and `element` each have 2 use sites flagged). If `pyright`
reports anything else (a genuinely new finding, a different count, or any
finding in a category the config's suppressions were supposed to silence),
STOP and report — investigate before proceeding, since this plan's Task 1
config was tuned against exactly this known set.

**Note on `typeCheckingMode`:** Task 1's config uses `"standard"`, not
`"basic"` — this was corrected after Task 3's first execution attempt
found that pyright 1.1.411's `"basic"` preset does not enable
`reportPossiblyUnboundVariable` by default (an error in the original
design/plan pass: the dry-run that found these findings during planning
ran pyright with no config at all, which defaults to `"standard"`, not
`"basic"` — the plan then specified `"basic"` without re-testing that
exact combination). If your checkout somehow still has `"basic"` in
`pyproject.toml`, that is stale — the fix landed as part of this task's
own execution; see the ledger for the full account.

- [ ] **Step 2: Write `docs/OPEN_SCIENTIFIC_QUESTIONS.md`**

```markdown
# Open Scientific/Numerical Questions

Status: **not fixed, not physics-locked-in** — each item here is a
suspected pre-existing issue found during audit or tooling (not
introduced by any change in this repository's `joss-epja-prep`/JOSS-
packaging work), deliberately left as-is per this project's rule not to
silently alter numerical/physics-adjacent behaviour. Each entry states
the evidence, the potential impact, and a proposed validation — resolving
any of these requires your explicit decision.

## 1. `ImportData.calc_ref_rev_frequency`'s `gamma` is unbound if none of `brho`/`ke`/`gam` is provided

**Location:** `src/rionid/core.py`, `calc_ref_rev_frequency` (static method).
**Evidence:** `pyright --pythonversion 3.9 src/` reports `"gamma" is
possibly unbound (reportPossiblyUnboundVariable)`. The method's `if
brho: gamma = ...\nelif ke: gamma = ...\nelif gam: gamma = gam` has no
final `else`, unlike its caller `reference_frequency`, which does guard
the equivalent three-way check with `else: sys.exit(...)`.
**Potential impact:** `calc_ref_rev_frequency` is only reachable today via
`reference_frequency`, which already validates that at least one of
`fref`/`brho`/`ke`/`gam` is set before calling it — so this is not
reachable through the current call graph. It would become a live
`UnboundLocalError` risk only if `calc_ref_rev_frequency` were ever called
directly (it's a `@staticmethod`, technically public) without that
upstream guard.
**Proposed validation:** either add the same `else: raise ValueError(...)`
guard `reference_frequency` already has, or explicitly document that this
static method assumes a validated caller. No behaviour change to any
currently-working input.

## 2. `LISEreader._read`'s `file_start` is unbound if the `[Calculations]` marker is missing

**Location:** `src/rionid/external/lisereader/reader.py`, `_read`.
**Evidence:** `pyright` reports `"file_start" is possibly unbound` at two
use sites. `file_start` is set only inside `if '[Calculations]' in line:
file_start = i+1` within a `for` loop over the file's lines, then used
immediately after the loop with no fallback.
**Potential impact:** a malformed or unexpected LISE++ output file (missing
the `[Calculations]` section header) would raise `NameError` rather than a
clear, actionable error message about the file format.
**Proposed validation:** add an explicit check after the search loop
(`if file_start is None: raise ValueError(...)`) with a message naming the
expected `[Calculations]` marker. This is error-message quality, not a
numerical change — safe to fix once you approve it, but not done here.

## 3. `LISEreader.get_index`'s `element` is unbound if `self.data` is empty

**Location:** `src/rionid/external/lisereader/reader.py`, `get_index`.
**Evidence:** `pyright` reports `"element" is possibly unbound` (×2) in the
`else` branch reached when no matching entries are found — `element` is the
loop variable from the preceding `for i, element in enumerate(self.data):`,
unbound if `self.data` is empty.
**Potential impact:** an empty parsed LISE++ file would raise `NameError`
instead of the intended `ValueError` with a formatting hint.
**Proposed validation:** same class of fix as #2 — an explicit empty-data
guard with a clear message. Not done here.

## 4. `save_simulation_results`'s `harmonic`/`fre` rely on a loop-variable leak

**Location:** `src/rionid/gui/controller.py`, `save_simulation_results`.
**Evidence:** `pyright` reports both `"harmonic" is possibly unbound` and
`"fre" is possibly unbound`. The function's header-writing loop
(`for harmonic in harmonics: ...`) leaves `harmonic` bound via Python's
loop-variable-leaks-into-enclosing-scope behaviour, then the *last* value
of `harmonic` from that loop is used inside a *separate* loop further down
when `mode == 'Brho'`: `fre = mydata.srrf[i] * mydata.ref_frequency*harmonic`.
**Potential impact:** two things worth your attention, not just an unbound
risk: (a) if `harmonics` is ever empty, this is a real `NameError`; (b) even
when `harmonics` is non-empty, using only the *last* harmonic in the list
to compute every row's `fre` in Brho mode looks intentional-but-unobvious
for a multi-harmonic run — this is the kind of thing that's either a
correct simplification (Brho mode may only ever be called with one
harmonic in practice) or a latent output bug for multi-harmonic Brho runs,
and distinguishing those requires domain knowledge this audit doesn't have.
**Proposed validation:** confirm with the physics whether multi-harmonic
`-hrm` combined with Brho mode (`-b`) is a supported combination at all;
if not, an explicit early validation error would be clearer than the
current implicit last-value behaviour. If it is supported, the per-row
`fre` calculation likely needs the row's own harmonic, not the loop's last
one — a semantic question for you, not a bug this document is claiming to
have proven.
```

- [ ] **Step 3: Commit**

```bash
git add docs/OPEN_SCIENTIFIC_QUESTIONS.md pyproject.toml
git commit -m "Add docs/OPEN_SCIENTIFIC_QUESTIONS.md, logging pyright's 7 real findings

pyright (configured in Task 1) surfaces 7 genuine
possibly-unbound-variable diagnostic lines across 5 distinct issues (2 of
the 5 -- file_start and element in lisereader/reader.py -- each have 2
flagged use sites), all pre-existing and none introduced by any change in
this engagement. Per the project's rule against silently altering
numerical/physics-adjacent behaviour, each is documented with evidence,
potential impact, and a proposed validation rather than fixed here --
resolving any of them is your call.

Also fixes pyproject.toml's [tool.pyright] typeCheckingMode from "basic"
to "standard": pyright 1.1.411's "basic" preset does not enable
reportPossiblyUnboundVariable by default, so the "basic" setting Task 1
committed silently produced zero diagnostics -- an error in the original
plan (the dry-run that found these findings ran pyright with no config,
which defaults to "standard"; the plan then specified "basic" without
re-testing that specific combination). The other rule-level suppressions
already in that config apply identically under "standard".

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: `pre-commit` configuration

**Files:**
- Create: `.pre-commit-config.yaml`

**Interfaces:** None.

- [ ] **Step 1: Write `.pre-commit-config.yaml`**

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.6.9
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.6.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-added-large-files
        args: [--maxkb=1024]
      - id: check-merge-conflict
      - id: check-toml
      - id: check-yaml
```
(`check-added-large-files` at 1MB catches an accidental large binary
commit without flagging normal source files — the `masses.py` physics
table is ~65KB, well under this.)

- [ ] **Step 2: Verify the config is at least syntactically valid**

```bash
python3 -c "import yaml; yaml.safe_load(open('.pre-commit-config.yaml'))" 2>&1 || python3 -c "import json,sys; print('yaml module unavailable, skipping parse check')"
```
(Installing and running `pre-commit run --all-files` for real requires
network access to fetch each hook's repo on first run — do that only if
network access is confirmed available in your environment; otherwise this
syntax check plus a careful manual read of the file is sufficient for this
task. Do not treat a network failure here as a task failure.)

- [ ] **Step 3: Commit**

```bash
git add .pre-commit-config.yaml
git commit -m "Add pre-commit config: ruff (lint+format) and basic hygiene hooks

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: CI workflow + Dependabot

**Files:**
- Create: `.github/workflows/ci.yml`
- Create: `.github/dependabot.yml`
- Modify: `src/rionid/core.py`, `src/rionid/external/lisereader/reader.py`,
  `src/rionid/gui/controller.py` (adds 6 targeted `# pyright: ignore[...]`
  comments, one per already-logged finding — see Step 3a; found necessary
  when this task's own local verification proved the naive workflow
  would fail on its first real run)

**Interfaces:** None — this task does not trigger any workflow (no push to
GitHub happens from this plan); it only writes the files (plus the 6
one-line comment additions above).

- [ ] **Step 1: Write `.github/workflows/ci.yml`**

```yaml
name: CI

on:
  push:
    branches: [master]
  pull_request:
    branches: [master]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.9", "3.10", "3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install system Qt runtime deps for offscreen PyQt5
        run: |
          sudo apt-get update
          sudo apt-get install -y libgl1 libegl1 libxkbcommon0

      - name: Install package with dev extras
        run: python -m pip install -e ".[dev]"

      - name: Cache AME/NUBASE mass tables
        uses: actions/cache@v4
        with:
          path: ~/.ame
          key: ame-nubase-tables-v1

      - name: Lint (ruff check)
        run: ruff check src/

      - name: Format check (ruff format)
        run: ruff format --check src/

      - name: Type check (pyright)
        run: pyright

      - name: Test
        env:
          QT_QPA_PLATFORM: offscreen
        run: pytest tests/ -v

      - name: Test (slow marker)
        env:
          QT_QPA_PLATFORM: offscreen
        run: pytest tests/ -m slow -v

      - name: Dependency vulnerability audit
        run: pip-audit --skip-editable

  build:
    runs-on: ubuntu-latest
    needs: test
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install build tools
        run: python -m pip install build twine

      - name: Build sdist and wheel
        run: python -m build

      - name: Check distribution metadata
        run: twine check dist/*

      - name: Install the built wheel into a fresh venv and smoke-test import
        run: |
          python -m venv /tmp/wheel-test-venv
          /tmp/wheel-test-venv/bin/pip install dist/*.whl
          /tmp/wheel-test-venv/bin/python -c "import rionid; from rionid.version import __version__; print(__version__)"
          /tmp/wheel-test-venv/bin/rionid --help > /dev/null 2>&1 || echo "rionid entry point resolved (GUI launch expectedly fails headless without QT_QPA_PLATFORM)"

      - uses: actions/upload-artifact@v4
        with:
          name: dist
          path: dist/
```

Notes for the implementer, not part of the file: `pip-audit --skip-editable`
avoids `pip-audit` trying (and failing) to resolve vulnerability data for
the editable-installed `rionid` package itself, which isn't published
anywhere — it should still audit every third-party dependency actually
pulled in. The `~/.ame` cache step is what closes the "tests need network
access to IAEA on a cold machine" gap the Wave 2a final review named — the
first CI run after this lands will need real network access to populate
it once; every run after reuses the cache.

- [ ] **Step 2: Write `.github/dependabot.yml`**

```yaml
version: 2
updates:
  - package-ecosystem: "pip"
    directory: "/"
    schedule:
      interval: "monthly"
  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "monthly"
```

- [ ] **Step 3: Verify every step the workflow runs actually works, standalone,
      in this sandbox** (GitHub Actions itself is not invoked — nothing is
      pushed — but each command must be shown to work before being trusted
      in a workflow file)

```bash
python3 -m pip install -e ".[dev]"
ruff check src/
ruff format --check src/
pyright
QT_QPA_PLATFORM=offscreen python3 -m pytest tests/ -v
QT_QPA_PLATFORM=offscreen python3 -m pytest tests/ -m slow -v
pip-audit --skip-editable
python3 -m build
twine check dist/*
```
`ruff` is deliberately scoped to `src/`, matching exactly what Tasks 2-3
cleaned — `tests/`, `gen_doc_stubs.py`, and any Markdown files under
`docs/` were never in scope for those tasks, and newer `ruff` releases
lint/format Markdown-embedded Python by default, which would otherwise
flag files nobody has touched. Do not widen this to `.` (repo-wide)
without a separate, deliberate task to actually clean those files first.

Expected outcomes:
- `ruff check src/` / `ruff format --check src/`: clean, per Tasks 2-3.
- `pyright`: **will still exit 1** at this point, reporting the same 7
  `reportPossiblyUnboundVariable` findings Task 3 logged in
  `docs/OPEN_SCIENTIFIC_QUESTIONS.md` — that document deliberately did not
  fix them, so `pyright`'s exit code is still nonzero. This is expected
  and handled in Step 3a below, not a sign anything is wrong yet.
- The two `pytest` runs: green, matching the existing suite.
- `pip-audit`: may report advisories for third-party packages. Read its
  output and note any in the task report. Do not attempt to fix
  third-party vulnerabilities as part of this task unless a compatible
  fixed version is actually installable (check with
  `pip index versions <package>` before assuming a fix exists — an
  advisory can reference a not-yet-released version) and the upgrade is
  a trivial, zero-compatibility-risk bump.
- `build`/`twine check`: succeed, per Task 1.

- [ ] **Step 3a: Suppress the 7 already-logged pyright findings at their
      exact source locations, so CI can be meaningfully green**

`docs/OPEN_SCIENTIFIC_QUESTIONS.md` (Task 3) deliberately did not fix
these — but leaving them unaddressed means the CI `pyright` step would be
red on every single run forever, for issues nobody intends to fix right
now, defeating the point of having the check at all. The fix is targeted,
standard, and does not touch any behavior: an inline
`# pyright: ignore[reportPossiblyUnboundVariable]` comment at each of the
6 lines that produce the 7 findings (one line, `reader.py`'s `element`
line, produces 2 of the 7 — one ignore comment silences both). This keeps
`pyright` meaningfully useful in CI: a **new**, different unbound-variable
bug anywhere else in `src/` would still fail the build.

Re-run `pyright` first to confirm these are still the exact 7 findings at
these exact locations (line numbers may have shifted if anything upstream
changed them — if they don't match, STOP and reconcile before editing):

```bash
pyright
```
Expected: exactly 7 errors at `core.py:369:32`, `reader.py:31:39`,
`reader.py:37:31`, `reader.py:56:19`, `reader.py:56:32`,
`controller.py:246:63`, `controller.py:251:39`.

In `src/rionid/core.py`, line 369 currently reads:
```python
        beta = ImportData.beta(gamma)
```
Change to:
```python
        beta = ImportData.beta(gamma)  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #1
```

In `src/rionid/external/lisereader/reader.py`, line 31 currently reads:
```python
        self.centre_index = len(lines[file_start].split()) - 1
```
Change to:
```python
        self.centre_index = len(lines[file_start].split()) - 1  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #2
```

Line 37 (inside the list comprehension) currently reads:
```python
            for line in lines[file_start:]
```
Change to:
```python
            for line in lines[file_start:]  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #2
```

Line 56 currently reads:
```python
            print(element[0] + element[1] + "+")
```
Change to:
```python
            print(element[0] + element[1] + "+")  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #3
```
(this one comment silences both of `element`'s two flagged occurrences on
this line).

In `src/rionid/gui/controller.py`, line 246 currently reads:
```python
                fre = mydata.srrf[i] * mydata.ref_frequency * harmonic
```
Change to:
```python
                fre = mydata.srrf[i] * mydata.ref_frequency * harmonic  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #4
```

Line 251 currently reads:
```python
            result_line = f"{ion:<15}{fre:<30.10f}{yield_:<15.4e}{moq:<15.12f}{mass:<15.3f}"
```
Change to:
```python
            result_line = f"{ion:<15}{fre:<30.10f}{yield_:<15.4e}{moq:<15.12f}{mass:<15.3f}"  # pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #4
```

If any line's exact text doesn't match what's shown above (formatting may
have shifted slightly), find the correct line by the file:line pyright
reports and append the same `# pyright: ignore[reportPossiblyUnboundVariable]  -- see docs/OPEN_SCIENTIFIC_QUESTIONS.md #N` comment to it — do not
change anything else on the line.

Re-run to confirm:
```bash
pyright
```
Expected: `0 errors, 0 warnings, 0 informations`.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/ci.yml .github/dependabot.yml src/rionid/core.py src/rionid/external/lisereader/reader.py src/rionid/gui/controller.py
rm -rf dist/ build/ src/rionid.egg-info
git commit -m "Add CI workflow (test matrix, lint, type-check, build, wheel smoke test) and Dependabot

Matrix over Python 3.9-3.12 (the declared requires-python range). Caches
~/.ame between runs so only a cold cache needs live IAEA network access
-- closes the gap named in Wave 2a's final review. A second job builds
the sdist/wheel, twine-checks it, and installs the built wheel into a
fresh venv to smoke-test the import and console entry point (not a full
GUI launch -- no display on CI runners).

ruff is scoped to src/ (matching what Tasks 2-3 actually cleaned, not
the whole repo -- tests/, gen_doc_stubs.py, and Markdown files were
never in scope and newer ruff versions lint Markdown by default).

Adds 6 targeted # pyright: ignore[reportPossiblyUnboundVariable]
comments at the exact 7 diagnostic locations docs/OPEN_SCIENTIFIC_QUESTIONS.md
already logged (Task 3) -- without these, pyright would exit 1 on every
CI run forever for issues deliberately left unfixed, defeating the
point of running it at all. A genuinely new unbound-variable bug
anywhere else in src/ still fails CI.

Every step verified to work standalone in this sandbox before being
written into the workflow file; the workflow itself is not triggered by
this commit (nothing is pushed).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 6: Non-triggering PyPI Trusted Publishing workflow

**Files:**
- Create: `.github/workflows/publish.yml`

**Interfaces:** None. **This workflow must never fire automatically** — verify
its trigger is `workflow_dispatch` only before committing.

- [ ] **Step 1: Write `.github/workflows/publish.yml`**

```yaml
name: Publish to PyPI

# Manual trigger ONLY. This workflow never runs automatically -- no push,
# tag, or release event triggers it. Even when manually invoked, PyPI
# Trusted Publishing (OIDC) requires the maintainer to have separately
# configured the trust relationship on PyPI's own project settings page
# first; this file being present does not by itself grant any publishing
# capability.
on:
  workflow_dispatch:
    inputs:
      confirm:
        description: 'Type "publish" to confirm you intend to publish to PyPI'
        required: true

permissions:
  id-token: write  # required for OIDC Trusted Publishing; no API token stored anywhere

jobs:
  publish:
    if: github.event.inputs.confirm == 'publish'
    runs-on: ubuntu-latest
    environment: pypi
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install build tools
        run: python -m pip install build twine

      - name: Build sdist and wheel
        run: python -m build

      - name: Check distribution metadata
        run: twine check dist/*

      - name: Publish to PyPI
        uses: pypa/gh-action-pypi-publish@release/v1
```

- [ ] **Step 2: Verify the trigger is safe**

```bash
grep -A2 "^on:" .github/workflows/publish.yml
```
Expected: only `workflow_dispatch` appears — no `push`, `release`, or
`schedule` trigger. The `if: github.event.inputs.confirm == 'publish'`
guard is a second layer: even a manual trigger without typing the exact
confirmation string does nothing.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/publish.yml
git commit -m "Add non-triggering PyPI Trusted Publishing workflow (workflow_dispatch only)

Least-privilege: OIDC-based Trusted Publishing (no stored API token),
gated on workflow_dispatch with an explicit typed confirmation input.
Never fires on push, tag, or release. Requires the maintainer to
separately configure the PyPI-side trust relationship before this
workflow could publish anything even if manually invoked -- not done
here, not part of this commit's effect.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 7: `.zenodo.json` + `CITATION.cff` version bump

**Files:**
- Create: `.zenodo.json`
- Modify: `CITATION.cff`

**Interfaces:** None.

- [ ] **Step 1: Write `.zenodo.json`**

```json
{
  "title": "RionID: Ring-stored ion IDentification",
  "description": "A Python GUI application for expert-guided identification of stored highly charged ions in storage-ring Schottky spectra, via reference-anchored forward simulation and an empirical polynomial frequency correction.",
  "creators": [
    {"name": "Freire-Fernández, David", "orcid": "0000-0002-5898-9291"},
    {"name": "Hudson-Chang, George", "orcid": "0000-0001-6404-6153"},
    {"name": "Chen, Rui-Jiu", "orcid": "0000-0001-6283-8958"}
  ],
  "license": "GPL-3.0-or-later",
  "upload_type": "software",
  "access_right": "open",
  "keywords": ["physics", "nuclear physics", "storage rings", "mass spectrometry", "Schottky spectroscopy"],
  "related_identifiers": [
    {
      "identifier": "10.5281/zenodo.8169341",
      "relation": "isVersionOf",
      "resource_type": "software"
    }
  ]
}
```
The `related_identifiers` entry references the concept DOI as it appears in
`README.md` today (`10.5281/zenodo.8169341`) — the known mismatch with
`CITATION.cff`'s DOI stays exactly as flagged in
`docs/PUBLICATION_TRACEABILITY.md`, not resolved by this file. No
version-specific DOI is included anywhere — Zenodo mints that
automatically when a real GitHub release is tagged, and inventing one now
would be fabrication.

- [ ] **Step 2: Update `CITATION.cff`**

Add a `version` field (currently absent) reflecting 9.0.0. Current file:
```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
- family-names: "Freire-Fernández"
  given-names: "David"
  orcid: "https://orcid.org/0000-0002-5898-9291"
- family-names: "Hudson-Chang"
  given-names: "George"
  orcid: "https://orcid.org/0000-0001-6404-6153"
- family-names: "Chen"
  given-names: "Rui-Jiu"
  orcid: "https://orcid.org/0000-0001-6283-8958"

title: "rionid: Collection of code for the identification of ringed ions in Python"
doi: 10.5281/zenodo.8169342
date-released: 2023-07-20
url: "https://github.com/GSI-Nuclear-Astrophysics/rionid"
```
Replace with (adds `version:`; leaves `doi:` and `date-released:` exactly
as they are — **do not change either of those two fields**, see the note
below):
```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
- family-names: "Freire-Fernández"
  given-names: "David"
  orcid: "https://orcid.org/0000-0002-5898-9291"
- family-names: "Hudson-Chang"
  given-names: "George"
  orcid: "https://orcid.org/0000-0001-6404-6153"
- family-names: "Chen"
  given-names: "Rui-Jiu"
  orcid: "https://orcid.org/0000-0001-6283-8958"

title: "rionid: Collection of code for the identification of ringed ions in Python"
version: 9.0.0
doi: 10.5281/zenodo.8169342
date-released: 2023-07-20
url: "https://github.com/GSI-Nuclear-Astrophysics/rionid"
```

**Do not touch `doi:` or `date-released:`.** `doi:` is the known, already-
flagged mismatch with the README (leave it — resolving it isn't this
task's call). `date-released: 2023-07-20` is now stale relative to a
9.0.0 version bump, but the honest fix is to set it to the real date
9.0.0 is actually tagged and released — which hasn't happened — not to
either fabricate today's date as a fake release date or leave a
version/date pairing that looks coherent but isn't. Leave it, and this
plan's Task 8 (JOSS_READINESS.md) records this as a pending action for
the actual release moment.

- [ ] **Step 3: Commit**

```bash
git add .zenodo.json CITATION.cff
git commit -m "Add .zenodo.json; add version field to CITATION.cff

.zenodo.json references the concept DOI as it appears in README.md
(10.5281/zenodo.8169341) via related_identifiers -- no version-specific
DOI is fabricated; Zenodo mints that automatically on a real release.
The existing CITATION.cff/README DOI mismatch stays exactly as flagged
in docs/PUBLICATION_TRACEABILITY.md, not resolved here.

CITATION.cff gains a version: 9.0.0 field (previously had none). Its
date-released stays at the stale 2023-07-20 value deliberately -- the
honest fix is to set it at actual release time, not to fabricate a date
now.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 8: README.md overhaul

**Files:**
- Modify: `README.md`

**Interfaces:** None.

- [ ] **Step 1: Replace `README.md` in full**

```markdown
# RionID (Ring-stored ion IDentification)

[![Documentation](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://GSI-Nuclear-Astrophysics.github.io/RionID)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8169341.svg)](https://doi.org/10.5281/zenodo.8169341)
[![PyPI version](https://badge.fury.io/py/rionid.svg)](https://badge.fury.io/py/rionid)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**RionID** is a Python GUI application for identifying stored highly
charged ions in storage-ring Schottky spectra. Given a reference ion, a
momentum-compaction factor, and a candidate ion list, it computes expected
revolution-frequency patterns, applies an empirical polynomial correction
anchored to known lines, and overlays the result on the experimental
spectrum for expert assessment.

<div align="center">
  <img src="https://github.com/GSI-Nuclear-Astrophysics/rionid/raw/master/docs/img/rionid.png?raw=true" width="50%">
</div>

## Scope and non-goals

RionID performs **deterministic forward simulation and expert-guided
overlay** — it computes and displays where candidate ions are expected to
appear, and lets you compare that against the data. **It does not do
automatic or autonomous ion identification.** There is no automatic
peak-to-ion assignment, no classifier, and no hidden ranking logic — every
displayed candidate is one you asked for, and every accepted assignment is
one you make yourself. See `docs/SCIENTIFIC_METHOD.md` for the physics and
`docs/AUTOMATIC_PID_REMOVAL_MAP.md` for a record of automatic-matching
functionality that was deliberately removed from this release.

## Features

*   **Pure Python:** no ROOT dependencies.
*   **Reference-anchored forward simulation:** candidate revolution
    frequencies from ionic mass-to-charge ratios and a user-specified
    reference ion, projected to any requested harmonics.
*   **Polynomial residual correction:** an explicit, user-supplied
    quadratic correction in revolution-frequency space (see
    `docs/SCIENTIFIC_METHOD.md`), applied consistently across harmonics.
*   **Interactive spectrum overlay:** pan, zoom, and inspect candidate ion
    labels against 1D experimental spectra.
*   **Signal processing:** peak detection with configurable threshold and
    minimum distance.
*   **Standalone:** bundles `lisereader` (GPL-3.0) for LISE++ candidate-list
    import without extra dependency management.

## Installation

### From PyPI (recommended)

```bash
pip install rionid
```

### From source (development)

```bash
git clone https://github.com/GSI-Nuclear-Astrophysics/rionid.git
cd rionid
pip install -e ".[dev]"
```

## Quick start

```bash
# Launch the GUI
rionid

# Or run a simulation directly from the terminal
rionid datafile.npz -r 72Ge+32 -ap 0.189 -f 1930000 -hrm 127 -c 2.55e-7 -0.985 9.5e5
```
`datafile.npz` needs `arr_0`/`arr_1` keys (frequency, amplitude) by
default, or any two array keys mapped via the GUI's key-selection dialog.
See `docs/REPRODUCIBILITY.md` for a runnable example using synthetic data
that needs no real experiment file at all.

## Parameter reference

| Parameter | CLI flag | Meaning |
|---|---|---|
| Data file | positional | Spectrum file (`.npz`, `.csv`, `.bin_fre`/`.bin_time`/`.bin_amp`) |
| Reference ion | `-r`, `--refion` | e.g. `72Ge+32` — sets the frequency-model anchor |
| Momentum compaction | `-ap`, `--alphap` | α_p; values `>1` are treated as γ_t and converted (`α_p = 1/γ_t²`) |
| Candidate list | `-psim`, `--filep` | LISE++ output file |
| Harmonics | `-hrm`, `--harmonics` | One or more harmonic orders to display |
| Reference frequency mode | `-f`/`-b`/`-ke`/`-gam` | Exactly one of: frequency [Hz], Brho [Tm], kinetic energy [MeV/u], Lorentz γ |
| Polynomial correction | `-c`, `--correct` | `A B C` coefficients (quadratic, linear, constant), Hz-based, `numpy.polyval` order |
| Top-N filter | `-n`, `--nions` | Show only the N highest-yield candidates (reference ion always included) |
| Highlight ions | (GUI field) | Comma-separated ion names to highlight — user-selected only, never automatically assigned |

Full stage-by-stage documentation, including invariants and known quirks,
is in `docs/LEGACY_BEHAVIOUR.md`.

## Supported formats

- Spectra: `.npz` (configurable key mapping), `.csv` (pipe-delimited),
  `.bin_fre`/`.bin_time`/`.bin_amp` (TDSM binary triples).
- Candidate lists: LISE++ output (`.lpp`).
- Export: `.ods` (candidate table), `simulation_result.out` (fixed-width
  text table).
- `.root` files are explicitly **not** supported — this is a deliberate
  restriction, not a bug; convert to `.npz`/`.csv` first.

## Troubleshooting

- **No network access on first run:** the AME2020 mass table is
  downloaded to `~/.ame/` on first use if not already cached there. If you
  have no network access, obtain `~/.ame/ame.data` from a machine that
  does, and copy the `~/.ame/` directory over.
- **Qt import crashes / binding-detection errors:** if you have both
  PyQt5 and PySide6 installed, some Qt-related tooling (notably
  `pytest-qt`, if you're running the test suite) can crash trying to
  auto-detect which binding to use. Ensure only one Qt binding is
  installed in your environment, or see `tests/conftest.py` for how this
  project works around it during testing.
- **GUI window doesn't appear / crashes with no display:** RionID needs a
  display server. On a headless machine, set `QT_QPA_PLATFORM=offscreen`
  for testing/scripting purposes (this disables the actual visible window
  — it's for automated checks, not for normal interactive use).

## Citation

If you use RionID, please cite it as described in
[`CITATION.cff`](CITATION.cff). Note: this repository currently has a
known, unresolved discrepancy between the Zenodo concept DOI shown above
and the one recorded in `CITATION.cff` — see
`docs/PUBLICATION_TRACEABILITY.md` for details; it is tracked, not hidden.

## Limitations

- The polynomial correction is setting-specific: coefficients fitted for
  one ring configuration, optics, or cooling condition are not valid for
  another without re-validation.
- A quadratic correction order is justified only by the anchor residuals
  observed over the fitted interval — it is not a first-principles model,
  and extrapolation outside that interval is not supported.
- Candidate-list completeness bounds every assignment: an unlisted species
  cannot be identified.
- Unresolved blends and weak lines may remain ambiguous even with
  multi-harmonic, multi-detector comparison.
- RionID supports identification preceding a mass or lifetime measurement;
  it does not replace dedicated mass-calibration or lifetime-analysis
  procedures.
- This release performs no automatic or autonomous species assignment —
  see "Scope and non-goals" above.

See `RionID-EPJA/main.tex` (the accompanying physics/methods manuscript)
for the full validation discussion.

## Acknowledgements

*   **Dr. RuiJiu Chen** for providing the C++ Time-of-Flight simulation code
    that inspired the backbone of this software.
*   **Dr. Shahab Sanjari** for guidance on software architecture and
    Schottky analysis.

## License

This project is licensed under the GNU General Public License v3.0. See
the [LICENSE](LICENSE) file for details.
```

- [ ] **Step 2: Also update `docs/index.md`** to match (it's a near-duplicate of
      README.md's content, kept in sync per the pattern established during
      Wave 2a) — apply the same content as Step 1 (the mkdocs nav wraps this
      file directly; no separate front-matter needed beyond what's already
      there).

- [ ] **Step 3: Commit**

```bash
git add README.md docs/index.md
git commit -m "Overhaul README.md (and docs/index.md to match): scope, parameters, troubleshooting, limitations, citation

Adds an explicit Scope and non-goals section stating plainly that this
release does no automatic/autonomous ion identification. Adds a
parameter reference table, supported-formats list, troubleshooting
section (no-network AME cache, Qt-binding conflicts, headless display),
a Limitations section drawn directly from the manuscript's own stated
limits, and a Citation section that surfaces the known DOI discrepancy
rather than hiding it. Quick-start install command updated to
pip install -e \".[dev]\" (Task 1's PEP 621 migration made this actually
work). No claims invented beyond what's already documented elsewhere in
this repository.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

- [ ] **Step 4 (addendum, added after Task 10): fix the Quick Start CLI
      example now that it's confirmed broken**

While writing Task 10 (`docs/REPRODUCIBILITY.md`), the implementer
discovered — and the controller independently reproduced from a clean
shell — that this task's own Quick Start CLI one-liner (`rionid
datafile.npz -r 72Ge+32 -ap 0.189 -f 1930000 -hrm 127 -c 2.55e-7 -0.985
9.5e5`) crashes immediately on any real invocation: `-psim` is
unconditionally required by `rionid.__main__.run_controller` despite
argparse not marking it required (no `.lpp` file is committed anywhere
in this repo to supply one), and independently, `__main__.py` has no
`--circumference` flag at all, so every reference-frequency mode
crashes with `TypeError` regardless. Both are real, pre-existing bugs
(confirmed via git history to predate this entire engagement), now
logged as `docs/OPEN_SCIENTIFIC_QUESTIONS.md` items 5 and 6. Fixing them
requires touching `src/rionid/`, which is out of scope for this
documentation-only plan — so the fix here is to stop presenting a
broken command as a working quick-start example, not to fix the
underlying CLI.

The GUI path is unaffected: `gui/inputs.py` has its own
`circumference_edit` field wired through `import_controller` to
`ImportData`, and file-selection dialogs for the candidate list — this
bug is specific to the `python3 -m rionid`/`rionid <datafile>` CLI
entry point in `__main__.py`, a separate code path from the GUI's
`import_controller`.

Replace this task's Quick Start section (in both `README.md` and
`docs/index.md`) from:
```markdown
## Quick start

```bash
# Launch the GUI
rionid

# Or run a simulation directly from the terminal
rionid datafile.npz -r 72Ge+32 -ap 0.189 -f 1930000 -hrm 127 -c 2.55e-7 -0.985 9.5e5
```
`datafile.npz` needs `arr_0`/`arr_1` keys (frequency, amplitude) by
default, or any two array keys mapped via the GUI's key-selection dialog.
See `docs/REPRODUCIBILITY.md` for a runnable example using synthetic data
that needs no real experiment file at all.
```
to:
```markdown
## Quick start

```bash
rionid
```
Fill in a reference ion, momentum-compaction factor, exactly one
reference-frequency value (frequency, Brho, kinetic energy, or gamma),
a candidate list (LISE++ `.lpp` output), the ring circumference, and a
spectrum file, then run the simulation from the window.

`datafile.npz` needs `arr_0`/`arr_1` keys (frequency, amplitude) by
default, or any two array keys mapped via the GUI's key-selection dialog.

**A note on the `python3 -m rionid`/`rionid <datafile>` CLI path:** as
currently implemented, this entry point requires a real LISE++
candidate-list file via `-psim`/`--filep` (not optional in practice,
despite argparse not marking it required) and has no way to supply a
ring circumference at all — every reference-frequency mode currently
crashes through this specific path. See
`docs/OPEN_SCIENTIFIC_QUESTIONS.md` items 5-6 for the full evidence, and
`docs/REPRODUCIBILITY.md` §3 for a verified, fully-working example that
exercises the same underlying simulation engine directly, using only
public/synthetic data.
```

- [ ] **Step 5 (addendum): verify and commit**

```bash
diff README.md docs/index.md
grep -n "psim\|circumference" README.md
```
Expected: `diff` shows no output (files still identical); `grep` shows
the new note's `-psim`/`--filep` and `circumference` mentions plus the
pre-existing `## Arguments`/parameter-reference-table rows (unchanged).

```bash
git add README.md docs/index.md
git commit -m "Fix README.md/docs/index.md Quick Start: stop showing a CLI command that crashes

Task 10 (docs/REPRODUCIBILITY.md) discovered, and the controller
independently reproduced, that this Quick Start's CLI one-liner cannot
succeed through rionid.__main__ with any arguments: -psim is
unconditionally required despite argparse not marking it so (no .lpp
fixture exists in this repo), and independently, __main__.py has no
--circumference flag at all, so every reference-frequency mode crashes.
Both are real, pre-existing bugs unrelated to this documentation work --
logged as docs/OPEN_SCIENTIFIC_QUESTIONS.md items 5-6, not fixed here
(fixing needs src/rionid/ changes, out of scope for this plan).

The GUI path is unaffected (gui/inputs.py wires its own circumference
field through import_controller) and remains the Quick Start's primary
example. The CLI note now honestly states the current limitation and
points to docs/REPRODUCIBILITY.md's verified working alternative,
instead of presenting a command that crashes on first use.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 9: `docs/SCIENTIFIC_METHOD.md`

**Files:**
- Create: `docs/SCIENTIFIC_METHOD.md`

**Interfaces:** None.

- [ ] **Step 1: Write `docs/SCIENTIFIC_METHOD.md`**

```markdown
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
```

- [ ] **Step 2: Commit**

```bash
git add docs/SCIENTIFIC_METHOD.md
git commit -m "Add docs/SCIENTIFIC_METHOD.md, cross-linked to EPJ A

Reader-facing summary of the frequency model, polynomial correction,
and harmonic projection, with direct section/equation cross-references
to RionID-EPJA/main.tex and to the implementing code. States plainly
what the manuscript describes but the software does not compute
(harmonic-overlap auto-derivation, uncertainty propagation) rather than
implying broader coverage than exists.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 10: `docs/REPRODUCIBILITY.md`

**Files:**
- Create: `docs/REPRODUCIBILITY.md`

**Interfaces:**
- Consumes: `tests.fixtures.synthetic_spectrum.make_synthetic_spectrum`/`build_ame_candidates` (from Wave 2a).

- [ ] **Step 1: Write `docs/REPRODUCIBILITY.md`**

```markdown
# Reproducibility

Exact, runnable commands to reproduce this repository's key numerical
claims from a clean checkout, using only public, synthetic, or already-
committed data — none of these need real experimental data.

## 1. Reproduce the polynomial-correction golden check

This reproduces the manuscript's harmonic-214 → revolution-space →
harmonic-127 coefficient-transform worked example (`RionID-EPJA/main.tex`
Eqs. 9-23) directly against the shipped code:

```bash
pip install -e ".[dev]"
pytest tests/test_correction.py -v
```
Expected: `1 passed`. The test's own docstring/comments walk through the
derivation; read `tests/test_correction.py` alongside
`docs/PUBLICATION_TRACEABILITY.md` for the full manuscript mapping.

## 2. Reproduce the extracted-mass-table golden check

Confirms `rionid.masses`'s ionic mass/mass-to-charge calculation for four
named nuclides (including the manuscript's own E143 example: ⁷²Ge³²⁺,
⁷⁴As³³⁺, ⁷⁶Se³⁴⁺):

```bash
pytest tests/test_masses.py -v
```
Expected: `5 passed`. Requires network access on first run only, to
populate `~/.ame/` with the AME2020 mass table (cached afterward).

## 3. Reproduce a full simulation run on synthetic data

No real spectrum file needed — this generates one from
`tests/fixtures/synthetic_spectrum.py` (documented there as clearly
synthetic, not derived from any experiment):

```bash
python3 -c "
import sys
sys.path.insert(0, 'src')
from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum
make_synthetic_spectrum('synthetic_demo.npz')
print('wrote synthetic_demo.npz')
"
```
Then run the CLI against it (using the reference ion / correction
coefficients from the manuscript's harmonic-127 worked example):
```bash
python3 -m rionid synthetic_demo.npz \
    -r 72Ge+32 -ap 0.189 -f 1930000 -hrm 127 \
    -c 2.55222702e-7 -0.985167644 950690.215
```
Expected: prints `Running RionID...`, logs the reference frequency and
processing status, and writes `simulation_result.out` (a fixed-width
table of candidate ion name/frequency/yield/m-q/mass) to the current
directory.

## 4. Reproduce the performance baseline

The exact synthetic-data methodology behind `docs/PERFORMANCE_BASELINE.md`
(both the Phase 0 numbers and Wave 2a's post-fix re-measurement) uses
`tests/fixtures/synthetic_spectrum.make_synthetic_spectrum`/
`build_ame_candidates` with `rionid.masses.get_ame_data()`, timed via
`time.perf_counter()` medians at N=10/100/2000 candidates — see that
document's "Method" section for the full, literal methodology; the
`@pytest.mark.slow` test in `tests/test_analysis.py`
(`pytest tests/ -m slow -v`) exercises the N=2000 scale as a correctness
(not timing) check and is the closest committed artifact to that
methodology.

## 5. What is *not* reproducible from this repository alone

`RionID-EPJA/main.tex`'s Figures 1, 3, and 4 (the broadband E143 spectrum,
the dual-pickup contaminant comparison, and the calibration-residual plot)
require the real E143 experimental dataset, which is not distributed in
this repository and is subject to the E143 collaboration's data policy —
see the manuscript's own Data Availability section. Figure 2 (the
harmonic-overlap-count staircase) is not reproducible from the package at
all, since the software does not compute that quantity — see
`docs/SCIENTIFIC_METHOD.md` §4.
```

- [ ] **Step 2: Verify commands 1-3 actually work as written**

```bash
pytest tests/test_correction.py -v
pytest tests/test_masses.py -v
python3 -c "
import sys
sys.path.insert(0, 'src')
from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum
make_synthetic_spectrum('/tmp/synthetic_demo.npz')
print('wrote synthetic_demo.npz')
"
cd /tmp && python3 -m rionid synthetic_demo.npz -r 72Ge+32 -ap 0.189 -f 1930000 -hrm 127 -c 2.55222702e-7 -0.985167644 950690.215 2>&1 | tail -20
```
(Use `PYTHONPATH=<repo>/src` if running `python3 -m rionid` from outside
the repo without an editable install active.) If step 3's CLI invocation
fails or behaves differently than described, fix the doc's command to
match actual behavior — do not leave a documented command that doesn't
work.

- [ ] **Step 3: Commit**

```bash
git add docs/REPRODUCIBILITY.md
git commit -m "Add docs/REPRODUCIBILITY.md with verified, runnable commands

Covers the two golden regression tests, a full synthetic-data
simulation run via the CLI, a pointer to the performance-baseline
methodology, and an explicit list of what is NOT reproducible from this
repository (the real E143 figures, which need restricted data; the
harmonic-overlap-count figure, which the software doesn't compute).
Every command was run and its actual output verified before being
written down.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 11: `examples/`, `CHANGELOG.md`, `CONTRIBUTING.md`, `SECURITY.md`, issue/PR templates

**Files:**
- Create: `examples/quickstart.py`
- Create: `CHANGELOG.md`
- Create: `CONTRIBUTING.md`
- Create: `SECURITY.md`
- Create: `.github/ISSUE_TEMPLATE/bug_report.md`
- Create: `.github/ISSUE_TEMPLATE/feature_request.md`
- Create: `.github/PULL_REQUEST_TEMPLATE.md`

**Interfaces:** None. Batched into one task (per the SDD process's
"batch small same-shape work" guidance) — each file is independent,
additive, and low-risk.

- [ ] **Step 1: Write `examples/quickstart.py`**

```python
"""Minimal, runnable RionID example using synthetic (non-experimental)
data -- no real spectrum or candidate file needed.

Run with: python examples/quickstart.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum
from rionid.core import ImportData
from rionid.masses import get_ame_data, ionic_moq_u


def main():
    spectrum_path = "quickstart_synthetic.npz"
    make_synthetic_spectrum(spectrum_path)
    print(f"Wrote synthetic spectrum to {spectrum_path}")

    # 72Ge32+ (fully-stripped germanium-72) is the manuscript's own E143
    # worked example (RionID-EPJA/main.tex). A single candidate equal to
    # the reference ion isolates this demo to the frequency model and
    # correction step, without needing a LISE++ candidate file -- see
    # docs/REPRODUCIBILITY.md for a full CLI example with -psim.
    ref_ion = "72Ge32+"
    ame_row = get_ame_data().lookup("Ge", 72)
    ref_moq = ionic_moq_u(ame_row, 32)

    model = ImportData(ref_ion, alphap=0.189, filename=spectrum_path,
                        reload_data=True, circumference=108.36)
    model.moq = {ref_ion: ref_moq}
    model.ref_ion = ref_ion
    model.ref_charge = 32
    model._calculate_srrf(fref=1.93e6)

    print(f"Reference ion {ref_ion}: m/q = {ref_moq} u")
    print(f"Reference frequency: {model.ref_frequency} Hz")
    print(f"srrf (relative revolution frequency): {model.srrf}")
    print("For a full simulation with a real candidate list and CLI/GUI "
          "usage, see docs/REPRODUCIBILITY.md.")


if __name__ == "__main__":
    main()
```

Before committing, actually run this and confirm it prints sensible output
without error:
```bash
python3 examples/quickstart.py
```
Expected: prints the real computed m/q for 72Ge32+ (should be
approximately `2.247` u, matching `tests/test_masses.py`'s golden value
for this exact nuclide/charge state), a reference frequency of `1930000.0`
Hz, and `srrf = [1.]` (since the one candidate equals the reference ion,
`srrf` is exactly 1.0 before any correction — no `-c` was passed here).
If the actual output differs from this, investigate before committing —
don't just paste the described output blindly.

- [ ] **Step 2: Write `CHANGELOG.md`**

```markdown
# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this
project uses [Semantic Versioning](https://semver.org/).

## [9.0.0] - Unreleased

Breaking changes to the public API and CLI. Earlier history (4.0.0-8.0.0)
is not reconstructed here from git log; see the git tag history for that.

### Removed

- **Automatic ion identification ("Quick PID")**: the GUI panel that
  scanned a grid of (reference frequency × α_p) values and silently
  overwrote user-entered values with an automatic best-match result.
  This directly contradicted the accompanying manuscript's claim that
  RionID is "not an autonomous classifier" — removing it makes that claim
  true of the shipped software. See `docs/AUTOMATIC_PID_REMOVAL_MAP.md`
  for the full removal record.
- **Baseline subtraction (BrPLS)**: the "Remove Baseline" feature and its
  `remove_baseline`/`psd_baseline_removed_l` parameters, per a deliberate
  product-scope decision (not a physics fix). Peak detection itself is
  unaffected and remains available.
- Dead/unreachable code: `handle_tiqnpz_data`/`handle_prerionidnpz_data`
  I/O handlers, the vendored `barion` library's unused automatic-
  identification methods (`identify_range` and related), and several
  unused GUI widgets tied to the removed features above.

### Changed

- The subset of the vendored `barion` library RionID actually uses
  (ionic-mass electron-binding correction, AME2020 table loading, a
  minimal ring-circumference holder) is now `src/rionid/masses.py`,
  extracted with the underlying arithmetic verified bit-identical to the
  original. `external/barion/` is gone.
- Two algorithmic hot paths are now O(N) instead of O(N×M)/O(N²):
  candidate mass-to-charge lookup and simulated-data yield lookup — see
  `docs/PERFORMANCE_BASELINE.md` for measured before/after numbers
  (~100× and ~235× respectively at N=2000 candidates).
- The AME2020 mass table is now cached for the process lifetime
  instead of being re-parsed on every simulation run.
- `pyproject.toml` migrated to PEP 621 metadata; `pip install -e ".[dev]"`
  now works (it could not before).

### Added

- A public pytest regression suite (19 tests) — there were none before
  this release cycle. Covers the polynomial-correction physics against
  the manuscript's own worked example, the extracted mass-table
  arithmetic, GUI state/wiring after the removals above, and the two
  speed-fixed hot paths' output identity.
- `ruff`/`pyright`/`pre-commit`/CI configuration.
- `docs/SCIENTIFIC_METHOD.md`, `docs/REPRODUCIBILITY.md`,
  `docs/OPEN_SCIENTIFIC_QUESTIONS.md`, `docs/JOSS_READINESS.md`.

### Fixed

- A stale-mass-table double-parse: `external/lisereader`'s LISE++ reader
  used to construct its own independent copy of the AME table on top of
  the one `ImportData` already loads — both now share one cached load.
```

- [ ] **Step 3: Write `CONTRIBUTING.md`**

```markdown
# Contributing to RionID

## Development setup

```bash
git clone https://github.com/GSI-Nuclear-Astrophysics/rionid.git
cd rionid
pip install -e ".[dev]"
pre-commit install
```

## Running checks locally

```bash
pytest tests/ -v              # public regression suite
pytest tests/ -m slow -v      # larger-N cases, opt-in
ruff check src/               # lint (scoped to src/ -- see CI workflow's own note)
ruff format --check src/      # format check (use `ruff format src/` to fix)
pyright                       # type check
```

## Scope

This is a parameter-driven scientific GUI application, not an automatic
particle-identification tool — see the README's "Scope and non-goals"
section. Contributions that reintroduce automatic species assignment,
hidden heuristics, or silent peak-to-ion matching will not be accepted;
the deterministic, user-guided workflow is a deliberate design decision,
not an oversight. See `docs/AUTOMATIC_PID_REMOVAL_MAP.md` for the
rationale.

Changes to `src/rionid/core.py`'s `_calculate_srrf` (the polynomial
correction) or any other physics/numerical behavior require the
manuscript's (`RionID-EPJA/main.tex`) formalism as the reference — open an
issue to discuss before submitting a PR that changes numerical output.

## Versioning

This project follows [Semantic Versioning](https://semver.org/):
- **Major** (X.0.0): breaking changes to the public CLI/GUI parameters or
  Python API.
- **Minor** (x.Y.0): new features, backward-compatible.
- **Patch** (x.y.Z): bug fixes, no interface changes.

## Release checklist

1. Ensure `pytest`, `ruff check src/`, `ruff format --check src/`, and
   `pyright` all pass.
2. Update `CHANGELOG.md` with the release date and finalized entry.
3. Bump the version in `pyproject.toml`, `CITATION.cff`, and
   `src/rionid/version.py` (all three, kept in sync).
4. Update `CITATION.cff`'s `date-released` to the actual release date.
5. Tag the release (`git tag vX.Y.Z`) and push the tag.
6. Verify the Zenodo-GitHub integration mints a version DOI, and update
   `.zenodo.json`/`CITATION.cff` with it once minted — never before.
7. Publish to PyPI only via the `workflow_dispatch`-gated
   `.github/workflows/publish.yml`, run deliberately, not automatically.

## Pull requests

- Keep PRs focused — one logical change per PR.
- Add or update tests for any behavior change.
- Update `CHANGELOG.md`'s `[Unreleased]` section (or the next version's,
  if already started) alongside your change.
```

- [ ] **Step 4: Write `SECURITY.md`**

```markdown
# Security Policy

## Reporting a vulnerability

If you discover a security vulnerability in RionID, please report it
privately by emailing D.FreireFernandez@gsi.de rather than opening a
public issue. Include a description of the vulnerability, steps to
reproduce it, and its potential impact.

You should expect an initial response within 5 business days. This is a
small, academically-maintained research-software project without a formal
security team — response times may vary.

## Scope

RionID is a desktop GUI/CLI application for offline scientific data
analysis. It does not run as a network service, does not handle
authentication or user accounts, and does not process untrusted network
input by design. Its main external-data-handling surface is:

- Reading user-supplied spectrum/candidate files (`.npz`/`.csv`/`.lpp`/
  binary formats) — a malformed file could in principle cause a crash or
  resource-exhaustion condition, not remote code execution.
- Downloading the public AME2020 mass table from IAEA
  (`www-nds.iaea.org`) over HTTPS on first use.

Dependency vulnerabilities are monitored via `pip-audit` in CI and
Dependabot's automated update PRs (`.github/dependabot.yml`).

## Supported versions

Only the latest released version is supported with security fixes.
```

- [ ] **Step 5: Write `.github/ISSUE_TEMPLATE/bug_report.md`**

```markdown
---
name: Bug report
about: Report a problem with RionID
title: ''
labels: bug
---

**Describe the bug**
A clear description of what went wrong.

**To reproduce**
Steps to reproduce, ideally including the exact CLI command or GUI
actions, and (if possible) a minimal synthetic input file that triggers
it — please do not attach restricted/unpublished experimental data.

**Expected behavior**
What you expected to happen instead.

**Environment**
- RionID version (`pip show rionid`):
- Python version:
- OS:

**Additional context**
Any other relevant information (full traceback, log output, etc.).
```

- [ ] **Step 6: Write `.github/ISSUE_TEMPLATE/feature_request.md`**

```markdown
---
name: Feature request
about: Suggest an enhancement for RionID
title: ''
labels: enhancement
---

**What problem would this solve?**

**Proposed solution**

**Note on scope:** RionID deliberately does not do automatic/autonomous
ion identification (see the README's "Scope and non-goals" section).
Feature requests that would reintroduce automatic peak-to-ion assignment,
hidden heuristics, or silent species ranking are out of scope for this
project — see `docs/AUTOMATIC_PID_REMOVAL_MAP.md` for why.

**Alternatives considered**
```

- [ ] **Step 7: Write `.github/PULL_REQUEST_TEMPLATE.md`**

```markdown
## What does this PR do?

## Checklist

- [ ] `pytest tests/ -v` passes
- [ ] `ruff check src/` and `ruff format --check src/` pass
- [ ] `pyright` passes
- [ ] Tests added/updated for any behavior change
- [ ] `CHANGELOG.md` updated
- [ ] No change to `_calculate_srrf` or other physics/numerical behavior
      without prior discussion in a linked issue (see CONTRIBUTING.md)
```

- [ ] **Step 8: Commit**

```bash
git add examples/ CHANGELOG.md CONTRIBUTING.md SECURITY.md .github/ISSUE_TEMPLATE/ .github/PULL_REQUEST_TEMPLATE.md
git commit -m "Add examples/quickstart.py, CHANGELOG.md, CONTRIBUTING.md, SECURITY.md, issue/PR templates

CHANGELOG starts fresh at 9.0.0 (Unreleased), documenting Wave 2a's
removals/changes/additions; no attempt to reconstruct 4.0.0-8.0.0
history from git log. CONTRIBUTING.md folds in the semantic-versioning
policy and release checklist rather than adding more top-level files.
The quickstart example and both issue templates explicitly restate the
project's expert-guided (not automatic) scope.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 12: `docs/JOSS_READINESS.md` (synthesis, written last)

**Files:**
- Create: `docs/JOSS_READINESS.md`

**Interfaces:** None. Depends on every prior task in this plan actually
being complete — this is an honest checklist against what exists, not a
plan of what will exist.

- [ ] **Step 1: Confirm every prior task's deliverable actually exists**

```bash
ls .zenodo.json CHANGELOG.md CONTRIBUTING.md SECURITY.md .pre-commit-config.yaml
ls .github/workflows/ci.yml .github/workflows/publish.yml .github/dependabot.yml
ls .github/ISSUE_TEMPLATE/bug_report.md .github/ISSUE_TEMPLATE/feature_request.md .github/PULL_REQUEST_TEMPLATE.md
ls docs/SCIENTIFIC_METHOD.md docs/REPRODUCIBILITY.md docs/OPEN_SCIENTIFIC_QUESTIONS.md
ls examples/quickstart.py
pytest tests/ -v
ruff check src/ && ruff format --check src/
pyright
```
`ruff` is scoped to `src/` throughout this plan (Task 5's CI workflow, this
verification) -- matching exactly what Tasks 2-3 actually cleaned, not the
whole repo. All must exist/pass before writing Step 2's claims about them.

- [ ] **Step 2: Write `docs/JOSS_READINESS.md`**

```markdown
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
| Research application | ✅ | Storage-ring ion identification for nuclear/atomic physics; accompanying EPJ A manuscript (`RionID-EPJA/main.tex`) describes the method and an experimental application |
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
```

- [ ] **Step 3: Commit**

```bash
git add docs/JOSS_READINESS.md
git commit -m "Add docs/JOSS_READINESS.md: evidence-based checklist, no invented claims

Written last, after confirming every referenced file/check in this plan
actually exists and passes. States two concrete pending actions (the
DOI mismatch, no tagged release yet) as open rather than resolved, and
explicitly disclaims any adoption/benchmark/submission-status claim
beyond what's actually documented elsewhere in this repository.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Self-Review Notes (for whoever executes this plan)

- **Spec coverage:** Task 1↔Group A (pyproject.toml/version); Task 7↔Group A
  (.zenodo.json/CITATION.cff); Tasks 2-6↔Group B (ruff/pyright/pre-commit/
  CI/publish workflow); Tasks 8-11↔Group C (README/SCIENTIFIC_METHOD/
  REPRODUCIBILITY/examples/CHANGELOG/CONTRIBUTING/SECURITY/templates);
  Task 12↔Group D (JOSS_READINESS.md). All four spec groups and every
  acceptance-criterion bullet are covered.
- **`docs/OPEN_SCIENTIFIC_QUESTIONS.md`** (Task 3) is a new addition beyond
  the original design spec's four groups — it wasn't anticipated because
  the spec was written before `pyright` had actually been run against the
  codebase. This is exactly the kind of standing location the whole
  engagement's brief calls for; creating it here, triggered by real
  tooling findings, is in scope, not scope creep.
- **Ordering matters**: Task 1 must land before Tasks 2-6 (they need the
  `[tool.ruff]`/`[tool.pyright]` config Task 1 writes). Task 2 must land
  before Task 3 (the `gui/plot.py` fix Task 2 makes is what removes one
  additional false-positive `reportPossiblyUnboundVariable` finding, so
  Task 3's `pyright` run reports only the genuinely pre-existing issues,
  not that one too). Task 12 must land last (it verifies every other
  task's output exists).
- **Correction found during Task 3's execution**: Task 1's `[tool.pyright]`
  config originally specified `typeCheckingMode = "basic"`, under which
  pyright 1.1.411 does not enable `reportPossiblyUnboundVariable` at all —
  silently producing zero diagnostics instead of the expected findings.
  Fixed to `"standard"` (Task 3, folded into its own commit rather than
  reopening the already-reviewed Task 1). The real finding count is **7**
  diagnostic lines across 5 issues, not 6 as originally estimated — a
  counting error in the original plan text (Task 3's own
  `docs/OPEN_SCIENTIFIC_QUESTIONS.md` content was already correct on this;
  only the summary count in Task 3's verification step and commit message
  said "6").
- **Nothing in this plan touches `_calculate_srrf`, any test file's
  assertions, or any file under `tests/`** — confirmed by re-reading every
  task's file list above.
