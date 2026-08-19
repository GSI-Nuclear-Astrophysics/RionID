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
