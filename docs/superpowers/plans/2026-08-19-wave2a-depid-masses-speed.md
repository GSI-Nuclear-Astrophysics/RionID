# Wave 2a: Regression Tests, Automatic-PID Removal, `masses.py` Extraction, Baseline-Removal Removal, Speed Fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a public pytest regression suite against the current code, then use it to safely remove the automatic-PID ("Quick PID") surface, remove the baseline-subtraction feature, extract the used subset of the vendored `barion` library into a new `masses.py` (dropping `external/barion/` and the unused `requests` dependency), and apply three evidence-backed speed fixes — all without changing the canonical polynomial correction's arithmetic or any other physics.

**Architecture:** No new module layout beyond one new file (`src/rionid/masses.py`) and one new test package (`tests/`). `core.py`, `gui/inputs.py`, `gui/app.py`, `gui/controller.py`, `io.py` shrink in place; `baseline.py` and `external/barion/` are deleted outright. The full `REFACTORING_PLAN.md` module split (`correction.py`/`calibration.py`/`analysis.py`/etc.) is explicitly out of scope — deferred to a later "Wave 2b" as pure moves with no logic change.

**Tech Stack:** Python 3.9 (repo supports 3.9–3.12; this environment runs 3.9.7 — avoid syntax newer than 3.9), pytest, pytest-qt, PyQt5 + pyqtgraph (headless via `QT_QPA_PLATFORM=offscreen`), numpy, fortranformat.

**Spec:** `docs/superpowers/specs/2026-08-19-wave2a-depid-masses-speed-design.md` (read this first — it has the full rationale; this plan does not repeat it).

## Global Constraints

- **Never change `ImportData._calculate_srrf`'s arithmetic, coefficient ordering `(A,B,C)=(quadratic,linear,constant)`, or units.** Every task that touches `core.py` must leave this method's body byte-for-byte identical.
- **No dead code left behind**: every deletion in this plan removes the calling code, the config/UI surface, and the docs claim together, in the same task.
- **Floating-point assertions use explicit `rtol`/`atol`** chosen per quantity (`pytest.approx(..., rel=1e-12)` for pure arithmetic identities; exact equality where output must be bit-identical after a pure refactor) — never a bare default tolerance.
- **`external/barion/` is only deleted once `tests/test_masses.py` passes against the new `masses.py`** (Task 4) — never delete it first and hope the port is correct.
- **This sandbox has no `Xvfb`**; all GUI tests run via `QT_QPA_PLATFORM=offscreen` (confirmed working in `docs/PERFORMANCE_BASELINE.md`). Set this in `tests/conftest.py` before any `PyQt5` import.
- **`poetry` is not installed in this environment.** Where the plan says "install", use `pip install <package>` directly against whatever Python environment `pytest`/`python3` resolve to in this repo checkout — do not block a step on `poetry` being present. `pyproject.toml` stays the declarative source of truth regardless.
- Every task ends with `git add` + `git commit` for exactly the files that task touches — do not batch multiple tasks into one commit.

---

### Task 1: Dev tooling — `pytest-qt`, drop unused `requests`, pytest config

**Files:**
- Modify: `pyproject.toml`

**Interfaces:**
- Produces: a working `pytest` + `pytest-qt` installation any later task's tests can rely on; `[tool.pytest.ini_options]` with `testpaths = ["tests"]`.

- [ ] **Step 1: Edit `pyproject.toml`**

Remove this line from `[tool.poetry.dependencies]` (confirmed unused anywhere in `src/` during the Phase 0 audit — `external/barion/amedata.py` uses the standard-library `urllib.request`, not `requests`):
```toml
requests = "^2.31.0"        # Required for Barion to download mass tables
```

Add `pytest-qt` to `[tool.poetry.group.dev.dependencies]`, so that section reads:
```toml
[tool.poetry.group.dev.dependencies]
pytest = "^8.0.0"
pytest-qt = "^4.2.0"
black = "^25.0.0"
```

Add a new section anywhere after `[tool.poetry.scripts]`:
```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "--strict-markers"
markers = [
    "slow: larger-N performance/regression cases, not run by default",
]
```

- [ ] **Step 2: Install `pytest-qt` locally for this checkout**

Run: `pip install pytest-qt`
(If `poetry` is available in your environment instead, `poetry lock && poetry install --with dev` is the equivalent — but do not block on this if `poetry` is absent, per Global Constraints.)

- [ ] **Step 3: Verify**

Run: `python3 -c "import pytest_qt; print('pytest-qt OK')"`
Expected: prints `pytest-qt OK`, no `ModuleNotFoundError`.

Run: `python3 -c "import ast; ast.parse(open('pyproject.toml').read())" 2>/dev/null; python3 -c "import tomllib" 2>/dev/null || python3 -c "import toml; toml.load(open('pyproject.toml'))" ` — this just needs to not raise; if `tomllib`/`toml` aren't available, instead just visually confirm the file has no duplicate `[tool.poetry.group.dev.dependencies]` sections and no stray `requests` line: `grep -n "requests\|pytest-qt" pyproject.toml`.
Expected: only the `pytest-qt = "^4.2.0"` line matches; no `requests` line.

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml
git commit -m "Add pytest-qt dev dependency, pytest config; drop unused requests dependency

requests was declared as a runtime dependency but never imported anywhere
in src/ -- external/barion/amedata.py uses urllib.request for its AME-table
download instead. Confirmed by grep during the Phase 0 audit
(docs/LEGACY_BEHAVIOUR.md).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: Synthetic test fixtures

**Files:**
- Create: `tests/__init__.py` (empty)
- Create: `tests/fixtures/__init__.py` (empty)
- Create: `tests/fixtures/synthetic_spectrum.py`
- Create: `tests/conftest.py`

**Interfaces:**
- Produces: `tests.fixtures.synthetic_spectrum.make_synthetic_spectrum(path, n_points=2000, seed=42) -> path`, `tests.fixtures.synthetic_spectrum.build_ame_candidates(ame_table, n, zz_max=92) -> list[tuple]` (tuple shape: `(element_name: str, aa: int, zz: int, nn: int, charge_states: list[int], yield_: float)`, the exact shape `ImportData._calculate_moqs`/`_simulated_data` expect in `particles_to_simulate`). Pytest fixtures `qapp` (session-scoped `QApplication`, offscreen) and `synthetic_spectrum_path` (function-scoped, writes a fresh temp `.npz`).

- [ ] **Step 1: Create `tests/__init__.py` and `tests/fixtures/__init__.py`**

Both empty files (makes `tests` and `tests.fixtures` importable packages, needed so `from tests.fixtures.synthetic_spectrum import ...` works from test modules).

- [ ] **Step 2: Write `tests/fixtures/synthetic_spectrum.py`**

```python
"""Synthetic (non-experimental) test fixtures for RionID's regression
suite. All data here is generated, not derived from any real measurement,
and is safe to commit/redistribute publicly.
"""
import numpy as np


def make_synthetic_spectrum(path, n_points=2000, seed=42):
    """Writes a small synthetic Schottky-like spectrum to `path` as a
    2-array .npz (arr_0=frequency [Hz], arr_1=amplitude [a.u.]), readable
    by rionid.io.handle_spectrumnpz_data's default keys.

    Not real experimental data: a flat noise floor plus three Gaussian
    lines on a 1 MHz synthetic band, used only to exercise the peak
    detection / plotting code paths deterministically and quickly in CI.
    """
    rng = np.random.default_rng(seed)
    freq = np.linspace(244.5e6, 245.5e6, n_points)
    amp = np.abs(1e-3 * (1 + 0.2 * rng.standard_normal(n_points)))
    for f0, height, width in [
        (244.71e6, 1.0, 800.0),
        (244.98e6, 0.6, 600.0),
        (245.22e6, 0.35, 700.0),
    ]:
        amp += height * np.exp(-0.5 * ((freq - f0) / width) ** 2)
    np.savez(path, arr_0=freq, arr_1=amp)
    return path


def build_ame_candidates(ame_table, n, zz_max=92):
    """Builds `n` candidate tuples in the exact shape
    `rionid.core.ImportData._calculate_moqs`/`_simulated_data` expect in
    `particles_to_simulate`: `(element_name, aa, zz, nn, [charge_states],
    yield)`. Built from REAL rows of a loaded AME table (real Z/N/A/element
    name; fully-stripped charge state) -- only "which N to pick" and the
    fixed yield=1.0 are synthetic, not the nuclide data itself.
    """
    candidates = []
    seen = set()
    for row in ame_table:
        nn, zz, aa, name = row[3], row[4], row[5], row[6]
        if zz < 1 or zz > zz_max:
            continue
        key = (name, aa)
        if key in seen:
            continue
        seen.add(key)
        candidates.append((name, aa, zz, nn, [zz], 1.0))
        if len(candidates) >= n:
            break
    return candidates
```

- [ ] **Step 3: Write `tests/conftest.py`**

```python
"""Shared pytest fixtures. QT_QPA_PLATFORM must be set to 'offscreen'
BEFORE any PyQt5 import happens anywhere in the test session -- this
module is pytest's first import, so it goes here.
"""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402


@pytest.fixture(scope="session")
def qapp():
    """A session-wide QApplication. Runs headlessly via the offscreen
    platform plugin -- confirmed working without Xvfb during the Phase 0
    audit (see docs/PERFORMANCE_BASELINE.md)."""
    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


@pytest.fixture
def synthetic_spectrum_path(tmp_path):
    from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum

    path = tmp_path / "synthetic_spectrum.npz"
    make_synthetic_spectrum(str(path))
    return str(path)
```

- [ ] **Step 4: Verify the fixtures work standalone**

Run:
```bash
python3 -c "
import sys; sys.path.insert(0, '.')
from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum, build_ame_candidates
import numpy as np
make_synthetic_spectrum('/tmp/rionid_test_spectrum.npz')
d = np.load('/tmp/rionid_test_spectrum.npz')
assert d['arr_0'].shape == (2000,) and d['arr_1'].shape == (2000,)
print('spectrum fixture OK')
"
```
Expected: prints `spectrum fixture OK`, no exception. (`build_ame_candidates` is exercised for real in Task 4 onward, once an `AMEData`-like table is available — no need to hand-construct a fake table here.)

- [ ] **Step 5: Commit**

```bash
git add tests/
git commit -m "Add synthetic test fixtures and shared pytest conftest

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: Golden regression test for the polynomial correction

**Files:**
- Create: `tests/test_correction.py`

**Interfaces:**
- Consumes: `rionid.core.ImportData` (unmodified in this task; its `_calculate_srrf`/`_simulated_data` methods, called directly via `ImportData.__new__(ImportData)` to avoid needing a full spectrum/candidate setup).
- Produces: a regression test that must keep passing, unmodified, through every subsequent task in this plan — it is the tripwire for any accidental change to the correction's arithmetic.

- [ ] **Step 1: Write `tests/test_correction.py`**

```python
"""Golden regression test for the canonical polynomial correction
(ImportData._calculate_srrf + ImportData._simulated_data), cross-checked
against RionID-EPJA/main.tex's own worked numeric example (Eqs. 9-21: the
harmonic-214 -> revolution-space -> harmonic-127 coefficient transform).

See docs/PUBLICATION_TRACEABILITY.md for the manuscript mapping. This test
must keep passing, UNMODIFIED, through every later task in this plan --
none of them may change _calculate_srrf's arithmetic, coefficient
ordering, or units.
"""
from types import SimpleNamespace

import pytest

from rionid.core import ImportData


def test_srrf_correction_matches_manuscript_harmonic_transform():
    # Manuscript's harmonic-214 coefficients (RionID-EPJA/main.tex:372-377)
    a214, b214, c214 = 1.19262945e-9, -9.85167644e-1, 2.03447706e8
    h214 = 214

    # Eq. (14)/(coefficient_h): revolution-space (A, B, C) from a
    # harmonic-h representation.
    A = h214 * a214
    B = b214
    C = c214 / h214
    assert A == pytest.approx(2.55222702e-07, rel=1e-8)
    assert B == pytest.approx(-9.85167644e-01, rel=1e-8)
    assert C == pytest.approx(9.50690215e05, rel=1e-8)

    # Eq. (coefficient_transfer): transform to harmonic 127.
    h127 = 127
    a127 = A / h127
    b127 = B
    c127 = C * h127
    assert a127 == pytest.approx(2.00962758e-09, rel=1e-8)
    assert c127 == pytest.approx(1.20737657e08, rel=1e-8)

    # Now exercise the REAL code path: ImportData._calculate_srrf
    # (Eq. rev_correction/corrected_rev) then ImportData._simulated_data
    # (Eq. harmonic, F=h*f), for a single candidate equal to the reference
    # ion. This isolates srrf to exactly f^(0)/f_ref BEFORE correction
    # (moq[candidate] - moq[ref] == 0), so f^(0) == f_ref exactly and the
    # test can pick f_ref freely without needing a real spectrum, candidate
    # list, or AME table load. See docs/superpowers/specs/
    # 2026-08-19-wave2a-depid-masses-speed-design.md for the derivation.
    model = ImportData.__new__(ImportData)
    model.moq = {"REF+1": 1.0}
    model.ref_ion = "REF+1"
    model.ref_charge = 1
    model.alphap = 0.189
    model.ring = SimpleNamespace(circumference=108.36)
    model.particles_to_simulate = []

    f0_test = 1.0e6  # Hz, arbitrary test revolution frequency
    model._calculate_srrf(fref=f0_test, correct=[A, B, C])

    expected_f_corr = f0_test + (A * f0_test**2 + B * f0_test + C)  # Eq. (corrected_rev)
    actual_f_corr = model.srrf[0] * model.ref_frequency
    assert actual_f_corr == pytest.approx(expected_f_corr, rel=1e-12)

    model._simulated_data(harmonics=[127.0], mode="Frequency")
    actual_F127 = float(model.simulated_data_dict["127.0"][0][0])

    # Direct harmonic-127 polynomial (Eq. h127_corrected), an independent
    # code path from the one above -- this is the actual scaling-invariant
    # check (Eq. scaling_invariant).
    F127_0 = f0_test * h127
    expected_F127 = F127_0 + a127 * F127_0**2 + b127 * F127_0 + c127
    # Matches to 0.000e+00 Hz per the Phase 0 audit; abs tolerance here is
    # generous float64 headroom, not a sign of expected disagreement.
    assert actual_F127 == pytest.approx(expected_F127, abs=1e-6)
```

- [ ] **Step 2: Run it**

Run: `python3 -m pytest tests/test_correction.py -v`
Expected: `1 passed`. (This is a characterization test capturing already-correct behavior, not red-green TDD — it must pass immediately against the current, unmodified `core.py`. If it fails, STOP and report the discrepancy rather than adjusting the test to match — this exact identity was independently verified during the Phase 0 audit.)

- [ ] **Step 3: Commit**

```bash
git add tests/test_correction.py
git commit -m "Add golden regression test for the polynomial correction

Ports the Phase 0 audit's manual verification (main.tex Eqs. 9-21,
harmonic-214 -> revolution-space -> harmonic-127 transform) into a
committed pytest test. Must keep passing unmodified through every later
task in this plan.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: Extract `masses.py` from vendored `barion`; delete `external/barion/`; O(1) AME lookup

**Files:**
- Create: `src/rionid/masses.py`
- Create: `tests/test_masses.py`
- Modify: `src/rionid/core.py`
- Modify: `src/rionid/external/lisereader/reader.py`
- Modify: `src/rionid/gui/controller.py` (imports `AMEData` from `external.barion` too — found during pre-flight review, not in the original design pass)
- Delete: `src/rionid/external/barion/` (whole directory: `amedata.py`, `particle.py`, `ring.py`, `__init__.py`, `__pycache__/`)

**Interfaces:**
- Produces: `rionid.masses.Ring(name, circumference)`, `rionid.masses.AMEData` (with `.ame_table`, `.nubase_table`, `.lookup(name, aa) -> row | None`, plus the `.to_mev`/`.to_u`/`.get_elbien` statics — `gui/controller.py`'s `save_simulation_results` calls `AMEData.to_mev` directly), `rionid.masses.get_ame_data() -> AMEData` (process-lifetime-cached singleton), `rionid.masses.ionic_mass_u(ame_row, qq) -> float`, `rionid.masses.ionic_moq_u(ame_row, qq) -> float`.
- Consumes (for the test only, before deletion): the still-present `rionid.external.barion.amedata.AMEData`/`rionid.external.barion.particle.Particle`, used once to derive the golden numbers baked into `tests/test_masses.py` below (already computed and verified during this planning pass — do not re-derive).

- [ ] **Step 1: Write `tests/test_masses.py` with pre-verified golden values**

These four nuclides' ionic mass/m-q were computed via the current (pre-extraction) `external.barion.particle.Particle` during planning, and are the ground truth this task's `masses.py` must reproduce exactly:

```python
"""Regression tests for src/rionid/masses.py -- the extracted, trimmed
subset of the vendored `barion` library RionID actually uses (ionic mass
with the electron-binding correction, AME table loading/lookup).

Golden nuclides: 12C6+ as a light, well-known cross-check, plus
72Ge32+/74As33+/76Se34+ because they are the manuscript's own worked E143
example (RionID-EPJA/main.tex). Expected values were computed via the
still-vendored external.barion.particle.Particle during the design of this
plan (docs/superpowers/plans/2026-08-19-wave2a-depid-masses-speed.md) --
do not re-derive them from masses.py itself, that would make this test
vacuous.
"""
import pytest

from rionid.masses import get_ame_data, ionic_mass_u, ionic_moq_u

# (element, mass_number, charge_state, expected_ionic_mass_u, expected_ionic_moq_u)
GOLDEN_NUCLIDES = [
    ("C", 12, 6, 11.996709621789908, 1.9994516036316512),
    ("Ge", 72, 32, 71.90458249215271, 2.2470182028797723),
    ("As", 74, 33, 73.90589141209702, 2.239572467033243),
    ("Se", 76, 34, 75.9006328732354, 2.232371555095159),
]


@pytest.fixture(scope="module")
def ame():
    return get_ame_data()


@pytest.mark.parametrize("name,aa,qq,expected_mass,expected_moq", GOLDEN_NUCLIDES)
def test_ionic_mass_and_moq(ame, name, aa, qq, expected_mass, expected_moq):
    row = ame.lookup(name, aa)
    assert row is not None, f"{aa}{name} not found in AME table"

    mass = ionic_mass_u(row, qq)
    moq = ionic_moq_u(row, qq)

    assert mass == pytest.approx(expected_mass, rel=1e-12)
    assert moq == pytest.approx(expected_moq, rel=1e-12)


def test_get_ame_data_is_cached():
    """The process-lifetime cache (docs/PERFORMANCE_BASELINE.md) must
    return the SAME object on repeated calls, not re-parse."""
    a = get_ame_data()
    b = get_ame_data()
    assert a is b
```

- [ ] **Step 2: Run it — must fail (masses.py doesn't exist yet)**

Run: `python3 -m pytest tests/test_masses.py -v`
Expected: `ModuleNotFoundError: No module named 'rionid.masses'` (or a collection error to that effect).

- [ ] **Step 3: Create `src/rionid/masses.py` — part 1 (head, via the Write tool)**

```python
"""Nuclide mass and mass-to-charge utilities.

Extracted from the vendored `barion` library (Xaratustrah, 2015-2016,
GPL-3.0 -- you are also an upstream co-owner of that library) to keep only
the subset RionID actually uses: AME2020/NUBASE2020 mass-table loading,
the ionic-mass electron-binding correction (RionID-EPJA/main.tex Sec. 2.1,
the f^(0) model, lines 198-200), and a minimal storage-ring circumference
holder. The physical constants and electron-binding-energy table
(`AMEData.ElBiEn`) below are copied verbatim from
`external/barion/amedata.py` and must not be hand-edited -- see
docs/AUTOMATIC_PID_REMOVAL_MAP.md (decisions #1-2) and REFACTORING_PLAN.md
for the extraction rationale.
"""
import os

import fortranformat as ff
import urllib.request as ur


class Ring:
    """A storage ring, reduced to what RionID uses: name and circumference.

    The upstream `barion.Ring` also carried `gamma_t`/`mag_rigidity`/
    `acceptance`/per-facility presets (`get_ring_dict`) -- confirmed unused
    by RionID in docs/LEGACY_BEHAVIOUR.md, dropped here.
    """

    def __init__(self, name, circumference):
        self.name = name
        self.circumference = circumference


class AMEData:
    """Loads and caches the AME2020 mass table and NUBASE2020 table.

    On construction, reads `~/.ame/ame.data` and `~/.ame/nubase.data`,
    downloading them first if absent (network side effect, documented in
    docs/LEGACY_BEHAVIOUR.md). Also builds an index for O(1) (element
    name, mass number) lookup via `lookup()`, replacing the linear table
    scan that was an O(N x table-size) hot spot -- see
    docs/PERFORMANCE_BASELINE.md.
    """

    AME_DATA_LINK = "https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt"
    AME_NUTAB_LINK = "https://www-nds.iaea.org/amdc/ame2020/nubase_3.mas20.txt"
    FOLDER_NAME = "/.ame/"

    def __init__(self):
        self.home_folder = os.path.expanduser("~") + AMEData.FOLDER_NAME
        os.makedirs(self.home_folder, exist_ok=True)
        self.ame_table = []
        self.nubase_table = []
        self.ame_data_filename = f"{self.home_folder}ame.data"
        self.nubase_data_filename = f"{self.home_folder}nubase.data"
        if not os.path.exists(self.ame_data_filename) or not os.path.exists(
            self.nubase_data_filename
        ):
            self._download()
        self._parse_ame()
        self._parse_nubase()
        self._index = {(row[6], row[5]): row for row in self.ame_table}

    def lookup(self, name, aa):
        """Returns the AME table row for (element name, mass number `aa`),
        or `None` if not present in the table."""
        return self._index.get((name, aa))

    def _download(self):
        """Downloads the AME2020/NUBASE2020 tables into `self.home_folder`.
        Ported verbatim from
        `external/barion/amedata.py:AMEData.download_ame_data`."""
        req = ur.Request(AMEData.AME_DATA_LINK, headers={"User-Agent": "Magic Browser"})
        g = ur.urlopen(req)
        with open(self.home_folder + "ame.data", "b+w") as f:
            f.write(g.read())

        req = ur.Request(AMEData.AME_NUTAB_LINK, headers={"User-Agent": "Magic Browser"})
        g = ur.urlopen(req)
        with open(self.home_folder + "nubase.data", "b+w") as f:
            f.write(g.read())

    def _parse_ame(self):
        """Parses the fixed-width AME2020 mass table. Ported verbatim from
        `external/barion/amedata.py:AMEData.init_ame_db`; the FORTRAN
        format string and 36-line header skip are the AME2020 file layout
        and must not change without a new AME release."""
        ffline = ff.FortranRecordReader(
            "(a4,a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6)"
        )
        with open(self.ame_data_filename) as f:
            for _ in range(36):
                next(f)
            for line in f:
                if "*" in line:
                    line = line.replace("*", "0")
                if "#" in line:
                    line = line.replace("#", ".")
                    line = "sys " + line
                else:
                    line = "exp " + line
                dataline = ffline.read(line)
                for i in range(len(dataline)):
                    if isinstance(dataline[i], str):
                        dataline[i] = dataline[i].strip()
                self.ame_table.append(dataline)

    def _parse_nubase(self):
        """Parses the fixed-width NUBASE2020 table. Ported verbatim from
        `external/barion/amedata.py:AMEData.init_nubase_db`."""
        with open(self.nubase_data_filename) as f:
            for _ in range(25):
                next(f)
            for line in f:
                name = line[11:16].strip()
                isomer = line[16:17].strip()
                if isomer != "":
                    continue
                lt = line[69:77].strip()
                mp = line[78:80].strip()
                self.nubase_table.append([name, lt, AMEData.get_multiplier(mp)])

    @staticmethod
    def to_mev(m_u):
        return m_u * AMEData.UU

    @staticmethod
    def to_u(m_mev):
        return m_mev / AMEData.UU

    @staticmethod
    def get_elbien(zz, qq):
        return AMEData.ElBiEn[zz][zz - qq]

    @staticmethod
    def get_multiplier(mult):
        """NUBASE half-life unit multiplier (seconds per unit). Ported
        verbatim from `external/barion/amedata.py:AMEData.get_multiplier`."""
        mult = mult.strip()
        table = {
            "": 0.0, "s": 1, "m": 60, "h": 3600, "d": 86400, "y": 31557600,
            "ms": 1e-3, "us": 1e-6, "ns": 1e-9, "ps": 1e-12, "fs": 1e-15,
            "as": 1e-18, "zs": 1e-21, "ys": 1e-24,
            "ky": 31557600e3, "My": 31557600e6, "Gy": 31557600e9,
            "Ty": 31557600e12, "Py": 31557600e15, "Ey": 31557600e18,
            "Zy": 31557600e21, "Yy": 31557600e24,
        }
        return table.get(mult, 0.0)
```

- [ ] **Step 4: Append the verbatim physics-data tail — mechanical copy, not hand-typed**

The electron-binding-energy table (`ElBiEn`, ~540 hardcoded numeric rows, cited to Rodrigues et al. 2004 / Huang et al. 1976 / Johnson & Soff 1988 / Plante et al. 1994) and the physical constants `CC`/`UU`/`EE`/`ME` must be copied byte-for-byte, never retyped. Run:

```bash
sed -n '262,837p' src/rionid/external/barion/amedata.py >> src/rionid/masses.py
```

This appends lines 262–837 of the still-present original file (the citation comment block, the four physical constants, and the full `ElBiEn = [...]` table) as a continuation of the `AMEData` class body — the appended block is already 4-space-indented at class-body level, so it must land directly after `get_multiplier`'s body with **no other code in between** (that's why this step runs immediately after Step 3, before any module-level code is added).

Verify the append landed correctly and matches the source exactly:
```bash
diff <(sed -n '262,837p' src/rionid/external/barion/amedata.py) <(tail -n 576 src/rionid/masses.py)
```
Expected: no output (files identical over that range). If this diff is non-empty, STOP — do not proceed with a mismatched physics table.

- [ ] **Step 5: Append the module-level helper functions**

These MUST be appended after Step 4's block (they are module-level, not class-body, so they cannot be interleaved with the class):

```bash
cat >> src/rionid/masses.py << 'PYEOF'


_ame_cache = None


def get_ame_data():
    """Returns a process-lifetime-cached `AMEData` instance.

    The AME/NUBASE tables are immutable for the life of the process (they
    only change if a user re-downloads them between runs), so re-parsing
    on every `ImportData` construction is wasted work -- see
    docs/PERFORMANCE_BASELINE.md, "AMEData() re-parses... on every
    ImportData construction".
    """
    global _ame_cache
    if _ame_cache is None:
        _ame_cache = AMEData()
    return _ame_cache


def ionic_mass_u(ame_row, qq):
    """Ionic mass in u for an ion in charge state `qq`, given its AME
    table row.

    Implements the manuscript's atomic-to-ionic mass correction
    (RionID-EPJA/main.tex:198-200, Sec. 2.1): the ionic mass is the atomic
    mass minus the removed electrons' rest mass, plus the corresponding
    change in total electron binding energy divided by c^2. Ported
    verbatim from
    `external/barion/particle.py:Particle.get_ionic_mass_in_u` -- do not
    change this arithmetic.
    """
    zz = ame_row[4]
    atomic_mass_u = ame_row[15] + ame_row[16] / 1.0e6
    if zz > 100:
        return atomic_mass_u
    return atomic_mass_u + AMEData.to_u(
        (AMEData.get_elbien(zz, 0) - AMEData.get_elbien(zz, qq)) / 1.0e6 - qq * AMEData.ME
    )


def ionic_moq_u(ame_row, qq):
    """Ionic mass-to-charge ratio in u, for an ion in charge state `qq`."""
    return ionic_mass_u(ame_row, qq) / qq
PYEOF
```

- [ ] **Step 6: Verify `masses.py` imports and parses cleanly**

Run: `python3 -c "import sys; sys.path.insert(0, 'src'); from rionid.masses import Ring, AMEData, get_ame_data, ionic_mass_u, ionic_moq_u; print('masses.py OK')"`
Expected: prints `masses.py OK`, no exception.

- [ ] **Step 7: Run `tests/test_masses.py` — must now pass**

Run: `python3 -m pytest tests/test_masses.py -v`
Expected: `5 passed` (4 parametrized nuclide cases + the cache-identity test).

- [ ] **Step 8: Update `src/rionid/core.py` to use `masses.py`, drop the dead `particles=` branch, and O(1)-index the AME lookup**

Replace the import block near the top of `core.py`:
```python
from rionid.external.barion.ring import Ring
from rionid.external.barion.amedata import AMEData
from rionid.external.barion.particle import Particle
from rionid.external.lisereader.reader import LISEreader
```
with:
```python
from rionid.masses import Ring, AMEData, get_ame_data, ionic_moq_u
from rionid.external.lisereader.reader import LISEreader
```

Replace `_set_particles_to_simulate_from_file`:
```python
    def _set_particles_to_simulate_from_file(self, particles_to_simulate):
        """Parses the LISE++ output file."""
        self.ame = AMEData()
        self.ame_data = self.ame.ame_table
        lise = LISEreader(particles_to_simulate)
        self.particles_to_simulate = lise.get_info_all()
```
with:
```python
    def _set_particles_to_simulate_from_file(self, particles_to_simulate):
        """Parses the LISE++ output file and loads the (process-cached)
        AME table -- see rionid.masses.get_ame_data."""
        self.ame = get_ame_data()
        self.ame_data = self.ame.ame_table
        lise = LISEreader(particles_to_simulate)
        self.particles_to_simulate = lise.get_info_all()
```

Replace `_calculate_moqs` in full:
```python
    def _calculate_moqs(self, particles = None):
        """Calculates mass-to-charge ratios for all particles."""
        self.moq = dict()
        self.total_mass = dict()
        
        if particles:
            for particle in particles:
                ion_name = f'{particle.tbl_aa}{particle.tbl_name}{particle.qq}+'
                m_q = particle.get_ionic_moq_in_u()
                self.moq[ion_name] = m_q
                self.total_mass[ion_name] = m_q * particle.qq
        else:
            for particle in self.particles_to_simulate:
                ion_name = f'{particle[1]}{particle[0]}{particle[4][-1]}+'
                for ame in self.ame_data:
                    if particle[0] == ame[6] and particle[1] == ame[5]:
                        pp = Particle(particle[2], particle[3], self.ame, self.ring)
                        pp.qq = particle[4][-1]
                        m_q = pp.get_ionic_moq_in_u()
                        self.moq[ion_name] = m_q
                        self.total_mass[ion_name] = m_q * pp.qq
                        self.protons[ion_name] = ame[4]
                        break
```
with:
```python
    def _calculate_moqs(self):
        """Calculates mass-to-charge ratios for every candidate in
        self.particles_to_simulate, via an O(1) AME-table lookup (see
        docs/PERFORMANCE_BASELINE.md; previously an
        O(N x AME-table-size) linear scan per candidate). The `particles=`
        parameter this method used to accept was dead code -- confirmed
        by repo-wide grep that every call site (__main__.py,
        gui/controller.py, gui/inputs.py) always calls it with no
        arguments -- and has been removed.
        """
        self.moq = dict()
        self.total_mass = dict()

        for particle in self.particles_to_simulate:
            ion_name = f'{particle[1]}{particle[0]}{particle[4][-1]}+'
            ame_row = self.ame.lookup(particle[0], particle[1])
            if ame_row is None:
                continue
            qq = particle[4][-1]
            m_q = ionic_moq_u(ame_row, qq)
            self.moq[ion_name] = m_q
            self.total_mass[ion_name] = m_q * qq
            self.protons[ion_name] = ame_row[4]
```

Every other `AMEData.CC` / `AMEData.UU` / `AMEData.to_mev` reference elsewhere in `core.py` (in `_calculate_srrf`, `calculate_brho_relativistic`, `gamma_brho`, `velocity`) needs no change — `AMEData` now resolves to `rionid.masses.AMEData`, which defines the same names with the same values (verbatim copy).

- [ ] **Step 9: Update `src/rionid/external/lisereader/reader.py`**

Replace:
```python
import numpy as np
from rionid.external.barion.amedata import AMEData
from re import sub

class LISEreader:
    def __init__(self, filename):
        ame = AMEData()
        ame.init_ame_db
        self.ame_data = ame.ame_table
        self._read(filename)
```
with:
```python
import numpy as np
from rionid.masses import get_ame_data
from re import sub

class LISEreader:
    def __init__(self, filename):
        # NOTE: previously constructed its own independent AMEData(),
        # re-parsing the mass table a SECOND time on top of core.py's own
        # load (docs/PERFORMANCE_BASELINE.md underestimated this -- the
        # real per-run AME-parse cost was ~2x the measured ~168ms). Now
        # shares the process-lifetime cache. The prior `ame.init_ame_db`
        # line (no parentheses) was a no-op that never actually called
        # the method -- dropped, since AMEData's own __init__ already
        # parses both tables.
        ame = get_ame_data()
        self.ame_data = ame.ame_table
        self._read(filename)
```

- [ ] **Step 9b: Update `src/rionid/gui/controller.py`**

`gui/controller.py` also imports `AMEData` from `external.barion` (used once, in `save_simulation_results`, to convert mass units for the output table) — found during the controller's pre-flight cross-task scan, not in the original design pass. Replace:
```python
from rionid.external.barion.amedata import AMEData
```
with:
```python
from rionid.masses import AMEData
```
No other change is needed in this file — `AMEData.to_mev(mass_u)` (line ~208) resolves identically against `masses.AMEData` (verbatim-copied constants/statics).

Verify no other reference to `external.barion` remains anywhere in `src/`:
```bash
grep -rln "external\.barion\|external import barion" src/
```
Expected after Steps 8, 9, and 9b: no matches.

- [ ] **Step 10: Delete `external/barion/`**

```bash
git rm -r src/rionid/external/barion/
rm -rf src/rionid/external/barion/__pycache__  # if git rm didn't already remove it (untracked cache dir)
```

- [ ] **Step 11: Run the full suite so far**

Run: `python3 -m pytest tests/ -v`
Expected: `tests/test_correction.py` (1 passed) + `tests/test_masses.py` (5 passed) all still green. If `test_correction.py` fails now, STOP — it means something in the import/refactor changed `_calculate_srrf`'s behavior, which must not happen.

- [ ] **Step 12: Update `README.md`'s "Standalone" bullet**

`barion` is no longer bundled as a separate vendored copy — only `lisereader` still is. Replace:
```markdown
*   **Standalone:** Bundles `barion` and `lisereader` (GPL-3.0) for easy installation without complex dependency management.
```
with:
```markdown
*   **Standalone:** Bundles `lisereader` (GPL-3.0) for easy installation without complex dependency management.
```

- [ ] **Step 13: Commit**

```bash
git add -A src/rionid/masses.py src/rionid/core.py src/rionid/gui/controller.py src/rionid/external/ tests/test_masses.py README.md
git commit -m "Extract used barion subset into masses.py; delete external/barion/; O(1) AME lookup

Extracts the ionic-mass electron-binding correction, AME/NUBASE table
loading, and a minimal ring-circumference holder from the vendored
barion library into src/rionid/masses.py -- the only subset RionID's own
code ever called (confirmed by call-graph audit in
docs/LEGACY_BEHAVIOUR.md). external/barion/ is deleted; its automatic
identify_range/get_unknown_*/get_nuclides_freqs methods and unused
enumeration helpers are not carried forward (per
docs/AUTOMATIC_PID_REMOVAL_MAP.md decisions #1-2, resolved with you as
an owner of the upstream barion library).

Also fixes _calculate_moqs's O(N x AME-table-size) linear scan (dict
lookup instead), drops its dead particles= parameter (never called with
an argument anywhere in the codebase), and de-duplicates a second,
previously-unnoticed AME-table parse that
external/lisereader/reader.py's LISEreader was independently performing
on every run -- both now share masses.get_ame_data()'s process-lifetime
cache. See docs/PERFORMANCE_BASELINE.md.

tests/test_masses.py's golden values were verified against the
still-vendored barion implementation before this change (see the design
spec) and now pass against masses.py directly.

gui/controller.py also imported AMEData from external.barion (used once
in save_simulation_results for a unit conversion) -- found during the
pre-flight cross-task review, repointed at masses.AMEData alongside the
rest.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Remove the automatic-PID ("Quick PID") surface

**Files:**
- Create: `tests/test_gui_smoke.py`
- Modify: `src/rionid/core.py`
- Modify: `src/rionid/gui/inputs.py`
- Modify: `src/rionid/gui/app.py`
- Modify: `README.md`
- Modify: `docs/index.md`

**Interfaces:**
- Consumes: `rionid.gui.inputs.RionID_GUI`, `rionid.gui.app.MainWindow` (both still constructible with no args / a `plot_widget=None` default, unaffected by this task's deletions).
- Produces: `RionID_GUI`/`MainWindow` with zero automatic-PID attributes, methods, or signals — the full "Remove" list of `docs/AUTOMATIC_PID_REMOVAL_MAP.md`.

- [ ] **Step 1: Write `tests/test_gui_smoke.py`**

```python
"""GUI regression tests: verify meaningful user-visible state and signal
wiring, not screenshots. Runs headlessly via QT_QPA_PLATFORM=offscreen
(tests/conftest.py sets this before any PyQt5 import).
"""


def test_quick_pid_surface_is_absent(qapp):
    """No automatic-PID controls, methods, or signals may exist on the
    input panel -- see docs/AUTOMATIC_PID_REMOVAL_MAP.md."""
    from rionid.gui.inputs import RionID_GUI

    panel = RionID_GUI()

    for attr in (
        "setup_quick_pid", "quick_pid_script", "onPlotClicked",
        "_stop_quick_pid", "overlay_sim_signal",
        "alphap_min_edit", "alphap_max_edit", "alphap_step_edit",
        "fref_min_edit", "fref_max_edit", "threshold_edit",
        "matched_result_edit",
    ):
        assert not hasattr(panel, attr), f"automatic-PID attribute still present: {attr}"


def test_collapsible_group_box_is_removed():
    import rionid.gui.inputs as inputs_module

    assert not hasattr(inputs_module, "CollapsibleGroupBox")


def test_compute_matches_is_removed():
    from rionid.core import ImportData

    assert not hasattr(ImportData, "compute_matches")
    assert not hasattr(ImportData, "save_matched_result")


def test_overlay_quick_pid_wiring_is_removed(qapp, qtbot):
    from rionid.gui.app import MainWindow

    win = MainWindow()
    qtbot.addWidget(win)
    assert not hasattr(win, "overlay_simulation")
    assert not hasattr(win.rion_input, "overlay_sim_signal")


def test_visualization_signal_still_updates_the_plot(qapp, qtbot):
    """The retained (non-automatic) run path's signal wiring must survive
    PID removal: RionID_GUI.visualization_signal ->
    MainWindow.update_visualization -> CreatePyGUI.updateData."""
    from rionid.gui.app import MainWindow

    win = MainWindow()
    qtbot.addWidget(win)

    calls = []
    win.visualization_widget.updateData = lambda data: calls.append(data)

    sentinel = object()
    win.rion_input.visualization_signal.emit(sentinel)

    assert calls == [sentinel]
```

- [ ] **Step 2: Run it — must fail (Quick PID still exists)**

Run: `python3 -m pytest tests/test_gui_smoke.py -v`
Expected: `test_quick_pid_surface_is_absent`, `test_collapsible_group_box_is_removed`, `test_compute_matches_is_removed`, `test_overlay_quick_pid_wiring_is_removed` all FAIL (assertion errors — the attributes/class still exist); `test_visualization_signal_still_updates_the_plot` PASSES already (that wiring isn't being touched).

- [ ] **Step 3: Remove `ImportData.compute_matches`/`chi2`/`match_count`/`save_matched_result` from `core.py`**

Delete the `compute_matches` method (`core.py`, originally lines 238–286) and the `save_matched_result` method (originally lines 288–295) in full.

In `ImportData.__init__`, remove these two lines from the "Results containers" block:
```python
        self.chi2 = 0
        self.match_count = 0
```
so that block reads just:
```python
        # Results containers
        self.peak_freqs = []
        self.peak_heights = []
```

- [ ] **Step 4: Remove Quick-PID UI and script from `gui/inputs.py`**

Delete the `setup_quick_pid` method in full, and its call inside `initUI`:
```python
        self.setup_file_selection()
        self.setup_parameters()
        self.setup_quick_pid()
        self.setup_controls()
```
becomes:
```python
        self.setup_file_selection()
        self.setup_parameters()
        self.setup_controls()
```

Delete the `threshold_edit` widget block from `setup_parameters` (it exists only to feed `quick_pid_script`):
```python
        # Threshold for PID
        self.threshold_edit = QLineEdit("1000")
        h_t = QHBoxLayout()
        h_t.addWidget(QLabel("Match Threshold (Hz):"))
        h_t.addWidget(self.threshold_edit)
        self.vbox.addLayout(h_t)
        
```
(delete this whole block; the blank line before `# Optional Features Group` stays).

Delete the `matched_result_edit` field from the "Optional Features" group in `setup_parameters`:
```python
        self.matched_result_edit = QLineEdit()
        opt_layout.addRow("Sim Result File:", self.simulation_result_edit)
        opt_layout.addRow("Matched Result File:", self.matched_result_edit)
```
becomes:
```python
        opt_layout.addRow("Sim Result File:", self.simulation_result_edit)
```
(keep `self.simulation_result_edit = QLineEdit()` above it — that field is retained, only `matched_result_edit` goes).

Delete `quick_pid_script` in full (the whole method).

Delete `onPlotClicked` in full:
```python
    @pyqtSlot()
    def onPlotClicked(self):
        self._stop_quick_pid = True
```

In `__init__`, remove:
```python
        self._stop_quick_pid = False
```

Remove the `overlay_sim_signal` and `clear_sim_signal` class-level signal declarations if `clear_sim_signal` is also unused — check first:
```bash
grep -n "clear_sim_signal" src/rionid/gui/*.py
```
If `clear_sim_signal` has no emitter/connector anywhere outside its declaration (expected — it was part of the same dead Quick-PID visual-feedback plumbing), remove both signal declarations:
```python
    visualization_signal = pyqtSignal(object)
    overlay_sim_signal = pyqtSignal(object)
    clear_sim_signal = pyqtSignal()
    signalError = pyqtSignal(str)
```
becomes:
```python
    visualization_signal = pyqtSignal(object)
    signalError = pyqtSignal(str)
```
(If `clear_sim_signal` DOES turn out to have another connector somewhere, keep it and only remove `overlay_sim_signal`.)

Delete the `CollapsibleGroupBox` class in full (unused, `gui/inputs.py`, originally lines 546–574 — this is the last class in the file).

Update `load_parameters` — remove the six Quick-PID-only lines:
```python
                self.matched_result_edit.setText(p.get('matched_result', ''))
                self.alphap_min_edit.setText(p.get('alphap_min', ''))
                self.alphap_max_edit.setText(p.get('alphap_max', ''))
                self.alphap_step_edit.setText(p.get('alphap_step', ''))
                self.fref_min_edit.setText(p.get('fref_min', ''))
                self.fref_max_edit.setText(p.get('fref_max', ''))
                self.threshold_edit.setText(p.get('threshold', '1000'))
```
(the line above them, `self.simulation_result_edit.setText(p.get('simulation_result', ''))`, stays — only these six go, and the `except FileNotFoundError: pass` line after stays).

Update `save_parameters` similarly — remove:
```python
            'matched_result': self.matched_result_edit.text(),
            'alphap_min': self.alphap_min_edit.text(),
            'alphap_max': self.alphap_max_edit.text(),
            'alphap_step': self.alphap_step_edit.text(),
            'fref_min': self.fref_min_edit.text(),
            'fref_max': self.fref_max_edit.text(),
            'threshold': self.threshold_edit.text()
```
so the dict literal's last two entries become:
```python
            'reload_data': self.reload_data_checkbox.isChecked(),
            'simulation_result': self.simulation_result_edit.text()
```
(remember to fix the trailing comma so this is still valid Python — `simulation_result` becomes the final entry, no trailing comma before the closing `}`).

**Do not touch** `_check_io_params` (the uncommitted NPZ-key-mapping helper) — leave it exactly as-is, even though `quick_pid_script` (its only current caller) is being deleted, per your earlier explicit instruction to leave that WIP untouched. It will simply be unused until you decide to wire it into `run_script` (which currently duplicates its logic inline) or remove it yourself later — flagged, not silently resolved, in this task's commit message.

- [ ] **Step 5: Remove the Quick-PID overlay wiring from `gui/app.py`**

Delete this line from `MainWindow.__init__`:
```python
        # 2. Overlay specific simulation (used for Quick PID visual feedback loop)
        self.rion_input.overlay_sim_signal.connect(self.overlay_simulation)
        
        # 3. Handle plot clicks (used to stop Quick PID loops or pick coordinates)
        self.visualization_widget.plotClicked.connect(self.rion_input.onPlotClicked)
```
so the signal-connections block becomes just:
```python
        # --- Signal Connections ---
        
        # Update plot when a full simulation run finishes
        self.rion_input.visualization_signal.connect(self.update_visualization)
```

Delete the `overlay_simulation` method in full:
```python
    def overlay_simulation(self, data):
        """
        Overlays a specific simulation result onto the existing plot.
        ...
        """
        self.visualization_widget.clear_simulated_data()
        self.visualization_widget.plot_simulated_data(data)
```

**Do not** remove `CreatePyGUI.plotClicked` in `gui/plot.py`, or `RionID_GUI.enterPlotPickMode`/`_onPlotPicked` in `gui/inputs.py` — these drive the retained manual "Pick" coordinate feature (`matching_freq_min_edit`/`matching_freq_max_edit`'s Pick buttons), per `docs/AUTOMATIC_PID_REMOVAL_MAP.md`'s Retain list.

- [ ] **Step 6: Update `README.md` and `docs/index.md`**

In both files, remove this bullet from "Features":
```markdown
*   **Automated Matching:** Includes Quick PID logic to scan $\alpha_p$ and Reference Frequency to find the best match ($\chi^2$ minimization).
```

- [ ] **Step 7: Run the GUI smoke test — must now pass**

Run: `python3 -m pytest tests/test_gui_smoke.py -v`
Expected: `5 passed`.

- [ ] **Step 8: Run the full suite**

Run: `python3 -m pytest tests/ -v`
Expected: all tests from Tasks 3–5 pass (`test_correction.py`, `test_masses.py`, `test_gui_smoke.py`).

- [ ] **Step 9: Commit**

```bash
git add src/rionid/core.py src/rionid/gui/inputs.py src/rionid/gui/app.py README.md docs/index.md tests/test_gui_smoke.py
git commit -m "Remove the automatic-PID (Quick PID) surface

Deletes ImportData.compute_matches/chi2/match_count/save_matched_result,
RionID_GUI.setup_quick_pid/quick_pid_script/CollapsibleGroupBox and their
widgets/config keys, and MainWindow's overlay_sim_signal/
overlay_simulation wiring -- the full 'Remove' list in
docs/AUTOMATIC_PID_REMOVAL_MAP.md. Removes the README/docs 'Automated
Matching' claim.

RionID_GUI._check_io_params (uncommitted WIP at the start of this
effort) is deliberately left untouched, even though quick_pid_script was
its only caller -- per earlier instruction to leave that WIP alone. It
is now unused; wiring it into run_script (which duplicates its logic
inline) or removing it is a follow-up decision, not made here.

CreatePyGUI.plotClicked and RionID_GUI.enterPlotPickMode/_onPlotPicked
(the manual coordinate-picking feature) are retained -- they are not
automatic PID.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 6: Remove the baseline-subtraction feature

**Files:**
- Delete: `src/rionid/baseline.py`
- Modify: `src/rionid/core.py`
- Modify: `src/rionid/gui/controller.py`
- Modify: `src/rionid/gui/inputs.py`
- Modify: `tests/test_gui_smoke.py`
- Modify: `README.md`

**Interfaces:**
- Produces: `ImportData.__init__` and `import_controller` with no `remove_baseline`/`psd_baseline_removed_l` parameters; `RionID_GUI` with no baseline-removal widgets.

- [ ] **Step 1: Extend `tests/test_gui_smoke.py`**

Append:
```python


def test_baseline_removal_surface_is_absent(qapp):
    from rionid.gui.inputs import RionID_GUI

    panel = RionID_GUI()
    for attr in ("remove_baseline_checkbox", "psd_baseline_removed_l_edit"):
        assert not hasattr(panel, attr), f"baseline-removal attribute still present: {attr}"


def test_baseline_module_is_removed():
    import importlib

    import pytest

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("rionid.baseline")
```

- [ ] **Step 2: Run it — must fail (baseline.py and its widgets still exist)**

Run: `python3 -m pytest tests/test_gui_smoke.py -v -k baseline`
Expected: both new tests FAIL.

- [ ] **Step 3: Delete `src/rionid/baseline.py`**

```bash
git rm src/rionid/baseline.py
```

- [ ] **Step 4: Update `src/rionid/core.py`**

Remove the import:
```python
from rionid.baseline import NONPARAMS_EST
```

Change `ImportData.__init__`'s signature — remove `remove_baseline`/`psd_baseline_removed_l`:
```python
    def __init__(self, refion, alphap, filename=None, reload_data=None, circumference=None,
                 highlight_ions=None, remove_baseline=False, psd_baseline_removed_l=1e6,
                 peak_threshold_pct=0.05, min_distance=10, matching_freq_min=None, 
                 matching_freq_max=None, io_params=None):
```
becomes:
```python
    def __init__(self, refion, alphap, filename=None, reload_data=None, circumference=None,
                 highlight_ions=None, peak_threshold_pct=0.05, min_distance=10,
                 matching_freq_min=None, matching_freq_max=None, io_params=None):
```

In the constructor body, remove:
```python
        self.remove_baseline = remove_baseline
        self.psd_baseline_removed_l = psd_baseline_removed_l
```
so that block reads:
```python
        self.peak_threshold_pct = float(peak_threshold_pct) if peak_threshold_pct else 0.05
        self.min_distance = float(min_distance) if min_distance else 10
        self.matching_freq_min = matching_freq_min
        self.matching_freq_max = matching_freq_max
        self.io_params = io_params or {} 
```

Replace the "NEW DATA PROCESSING BLOCK" (removing step 1, baseline removal, only):
```python
            # --- NEW DATA PROCESSING BLOCK ---
            if self.experimental_data is not None:
                freq, amp = self.experimental_data
                
                # 1. Baseline Removal
                if remove_baseline:
                    try:
                        est = NONPARAMS_EST(amp)
                        baseline = est.pls('BrPLS', l=psd_baseline_removed_l, ratio=1e-6)
                        amp = amp - baseline
                    except Exception as e:
                        print(f"Baseline removal failed: {e}")
                        traceback.print_exc()

                # 2. Log-Safety (Clip negatives)
                # Ensure all values are > 0 for logarithmic plotting. 
                # We use 1e-9 as a "floor" value.
                amp = np.maximum(amp, 1e-29)

                # 3. Normalization
                # Scale so the highest peak is 1.0
                max_val = np.max(amp)
                if max_val > 0:
                    amp = amp / max_val
                
                # Update the stored data
                self.experimental_data = (freq, amp)
            # ---------------------------------
```
with:
```python
            # --- Data processing: log-safety clip + normalise ---
            if self.experimental_data is not None:
                freq, amp = self.experimental_data

                # Log-Safety (Clip negatives): ensure all values are > 0
                # for logarithmic plotting. NOTE: the floor here (1e-29)
                # intentionally does not match gui/plot.py's separate
                # 1e-9 floor for the same purpose -- pre-existing,
                # unchanged, logged in docs/LEGACY_BEHAVIOUR.md, not
                # fixed as part of this removal.
                amp = np.maximum(amp, 1e-29)

                # Normalization: scale so the highest peak is 1.0
                max_val = np.max(amp)
                if max_val > 0:
                    amp = amp / max_val

                self.experimental_data = (freq, amp)
```

Remove the baseline-removal block from `_get_experimental_data` entirely:
```python
        # Baseline removal
        if self.remove_baseline and self.experimental_data:
            try:
                freq, psd = self.experimental_data
                est = NONPARAMS_EST(psd)
                baseline = est.pls('BrPLS', l=self.psd_baseline_removed_l, ratio=1e-6)
                self.experimental_data = (freq, psd - baseline)
            except Exception as e:
                traceback.print_exc()
```

Check whether `traceback` is still used anywhere else in `core.py`:
```bash
grep -n "traceback" src/rionid/core.py
```
If the only remaining hit is the `import traceback` line itself, remove that import too.

- [ ] **Step 5: Update `src/rionid/gui/controller.py`**

Remove `remove_baseline`/`psd_baseline_removed_l` from `import_controller`'s signature:
```python
def import_controller(datafile=None, filep=None, alphap=None, refion=None, harmonics=None, 
                      nions=None, amplitude=None, circumference=None, mode=None, value=None, 
                      reload_data=None, remove_baseline=False, psd_baseline_removed_l=1e6,
                      peak_threshold_pct=0.05, min_distance=10, highlight_ions=None, 
                      io_params=None, sim_scalingfactor=None, matching_freq_min=None, 
                      matching_freq_max=None, correct=None):
```
becomes:
```python
def import_controller(datafile=None, filep=None, alphap=None, refion=None, harmonics=None, 
                      nions=None, amplitude=None, circumference=None, mode=None, value=None, 
                      reload_data=None, peak_threshold_pct=0.05, min_distance=10, highlight_ions=None, 
                      io_params=None, sim_scalingfactor=None, matching_freq_min=None, 
                      matching_freq_max=None, correct=None):
```

Remove the two matching docstring parameter blocks (`remove_baseline : bool, optional` ... and `psd_baseline_removed_l : float, optional` ...).

In the `ImportData(...)` call, remove:
```python
                            remove_baseline=remove_baseline, psd_baseline_removed_l=psd_baseline_removed_l,
```
so the call becomes:
```python
        mydata = ImportData(refion, float(alphap), filename=datafile, reload_data=reload_data, 
                            circumference=circumference, highlight_ions=highlight_ions,
                            peak_threshold_pct=peak_threshold_pct, min_distance=min_distance,
                            matching_freq_min=matching_freq_min, matching_freq_max=matching_freq_max,
                            io_params=io_params)
```

- [ ] **Step 6: Update `src/rionid/gui/inputs.py`**

Remove the baseline widget block from `setup_parameters`:
```python
        # Baseline
        self.remove_baseline_checkbox = QCheckBox('Remove Baseline')
        self.vbox.addWidget(self.remove_baseline_checkbox)
        
        self.psd_baseline_removed_l_edit = QLineEdit("1000000")
        hb_bl = QHBoxLayout()
        hb_bl.addWidget(QLabel("Baseline l:"))
        hb_bl.addWidget(self.psd_baseline_removed_l_edit)
        self.vbox.addLayout(hb_bl)

```
(delete this whole block; the `# Alpha P` block that follows stays).

In `load_parameters`, remove:
```python
                self.remove_baseline_checkbox.setChecked(p.get('remove_baseline_checkbox', False))
                self.psd_baseline_removed_l_edit.setText(str(p.get('psd_baseline_removed_l', '1000000')))
```

In `save_parameters`, remove:
```python
            'remove_baseline_checkbox': self.remove_baseline_checkbox.isChecked(),
            'psd_baseline_removed_l': self.psd_baseline_removed_l_edit.text(),
```

In `run_script`, remove:
```python
        psd_l = self._get_float(self.psd_baseline_removed_l_edit, 1000000.0)
```
and, in the `args = argparse.Namespace(...)` call, remove:
```python
            remove_baseline=self.remove_baseline_checkbox.isChecked(),
            psd_baseline_removed_l=psd_l,
```
Since `psd_l` is now undefined, also update the `Namespace(...)` call's use of it — there is no other reference to `psd_l`, so simply deleting both lines above is sufficient; verify with:
```bash
grep -n "psd_l\b" src/rionid/gui/inputs.py
```
Expected after the edit: no matches.

- [ ] **Step 7: Update `README.md`**

Replace:
```markdown
*   **Signal Processing:** Built-in baseline subtraction (BrPLS) and peak detection.
```
with:
```markdown
*   **Signal Processing:** Peak detection with configurable threshold and minimum distance.
```

Remove these two lines from the "Arguments" list (confirmed during Task 4/5 planning that neither is actually a real CLI flag — `__main__.py`'s argparse never defined them, and the CLI path never threads `peak_threshold_pct` through to `ImportData` either; this was a pre-existing README/CLI mismatch, not something this removal creates):
```markdown
*   `--remove_baseline`: Apply baseline subtraction.
*   `--peak_threshold_pct`: Peak detection threshold (0.0 - 1.0).
```

- [ ] **Step 8: Run the GUI smoke test — must now pass**

Run: `python3 -m pytest tests/test_gui_smoke.py -v`
Expected: `7 passed`.

- [ ] **Step 9: Run the full suite**

Run: `python3 -m pytest tests/ -v`
Expected: all green.

- [ ] **Step 10: Commit**

```bash
git add -A src/rionid/baseline.py src/rionid/core.py src/rionid/gui/controller.py src/rionid/gui/inputs.py tests/test_gui_smoke.py README.md
git commit -m "Remove the baseline-subtraction feature (BrPLS)

Per your request: deletes src/rionid/baseline.py and every call site
(ImportData.__init__'s two baseline blocks, import_controller's
parameters, the GUI's 'Remove Baseline' checkbox and lambda field).
No pyproject.toml dependency changes -- baseline.py only used scipy,
which stays required for scipy.signal.find_peaks (peak detection is
unaffected and retained).

Also fixes a pre-existing README/CLI mismatch: --remove_baseline and
--peak_threshold_pct were documented as CLI flags but never defined in
__main__.py's argparse, and the CLI path never threaded
peak_threshold_pct through to ImportData either -- baseline removal was
already GUI-only. Peak detection itself (find_peaks, threshold, min
distance) is unaffected and remains GUI-configurable.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 7: Remove dead I/O handlers

**Files:**
- Modify: `src/rionid/io.py`
- Modify: `src/rionid/core.py`

**Interfaces:**
- No public interface change — `handle_spectrumnpz_data` (the only `.npz` handler actually reachable, per `docs/LEGACY_BEHAVIOUR.md`) is untouched.

- [ ] **Step 1: Remove `handle_tiqnpz_data` and `handle_prerionidnpz_data` from `io.py`**

Delete both functions in full (originally `io.py:72-101` and `io.py:124-141`).

- [ ] **Step 2: Update the import in `core.py`**

Replace:
```python
from rionid.io import (
    read_psdata, 
    handle_read_tdsm_bin, 
    handle_spectrumnpz_data, # probably I will just keep this option, delete the others
    handle_tiqnpz_data
)
```
with:
```python
from rionid.io import (
    read_psdata,
    handle_read_tdsm_bin,
    handle_spectrumnpz_data,
)
```

- [ ] **Step 3: Remove the commented-out dead TIQ dispatch in `_get_experimental_data`**

Replace:
```python
        elif ext == '.npz':
            self.experimental_data = handle_spectrumnpz_data(filename, **self.io_params)
            #if 'spectrum' in base:
            #    self.experimental_data = handle_spectrumnpz_data(filename, **self.io_params)
            #else:
            #    self.experimental_data = handle_tiqnpz_data(filename, **self.io_params)
```
with:
```python
        elif ext == '.npz':
            self.experimental_data = handle_spectrumnpz_data(filename, **self.io_params)
```

- [ ] **Step 4: Verify nothing else references the deleted functions**

Run: `grep -rn "handle_tiqnpz_data\|handle_prerionidnpz_data" src/ tests/`
Expected: no matches.

- [ ] **Step 5: Run the full suite**

Run: `python3 -m pytest tests/ -v`
Expected: all green (this task doesn't change any tested behavior — `handle_spectrumnpz_data` was already the only reachable `.npz` path).

- [ ] **Step 6: Commit**

```bash
git add src/rionid/io.py src/rionid/core.py
git commit -m "Remove dead handle_tiqnpz_data / handle_prerionidnpz_data I/O handlers

Confirmed unreachable: the .npz dispatch in core.py always called
handle_spectrumnpz_data; the TIQ branch was commented out and
handle_prerionidnpz_data was never imported. Per your decision
(docs/AUTOMATIC_PID_REMOVAL_MAP.md decision #3).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 8: Speed fix — `_simulated_data`'s O(N²) yield lookup

**Files:**
- Create: `tests/test_analysis.py`
- Modify: `src/rionid/core.py`

**Interfaces:**
- Consumes: `tests.fixtures.synthetic_spectrum.build_ame_candidates`, `rionid.masses.get_ame_data`.
- Produces: `ImportData._simulated_data` with identical output, O(N) instead of O(N²).

- [ ] **Step 1: Write `tests/test_analysis.py`**

```python
"""Golden-output regression tests for the O(N)/O(N^2) hot paths identified
in docs/PERFORMANCE_BASELINE.md. These tests capture output BEFORE the
speed fix in this task lands and must still pass, UNMODIFIED, after --
proving the fix is output-preserving, not just faster.
"""
import numpy as np
import pytest

from rionid.core import ImportData
from rionid.masses import get_ame_data
from tests.fixtures.synthetic_spectrum import build_ame_candidates


@pytest.fixture(scope="module")
def ame_table():
    return get_ame_data().ame_table


def _make_model(ame_table, synthetic_spectrum_path, n):
    candidates = build_ame_candidates(ame_table, n)
    ref_name, ref_aa, ref_zz = candidates[0][0], candidates[0][1], candidates[0][2]
    model = ImportData(
        f"{ref_aa}{ref_name}{ref_zz}+", alphap=0.189,
        filename=synthetic_spectrum_path, reload_data=True,
        circumference=108.36,
    )
    model.ame = get_ame_data()
    model.ame_data = model.ame.ame_table
    model.particles_to_simulate = candidates
    return model


@pytest.mark.parametrize("n", [10, 100])
def test_calculate_moqs_output(ame_table, synthetic_spectrum_path, n):
    model = _make_model(ame_table, synthetic_spectrum_path, n)
    model._calculate_moqs()

    assert len(model.moq) == n
    assert set(model.moq.keys()) == set(model.total_mass.keys())
    assert all(v > 0 for v in model.moq.values())


@pytest.mark.parametrize("n", [10, 100])
def test_simulated_data_yield_lookup_output(ame_table, synthetic_spectrum_path, n):
    model = _make_model(ame_table, synthetic_spectrum_path, n)
    model._calculate_moqs()
    model._calculate_srrf(fref=1.93e6)

    model._simulated_data(harmonics=[127.0], mode="Frequency")

    assert len(model.yield_data) == n
    # Every synthetic candidate carries yield=1.0 (build_ame_candidates).
    assert np.all(np.asarray(model.yield_data) == 1.0)
    assert model.simulated_data_dict["127.0"].shape == (n, 3)


@pytest.mark.slow
def test_simulated_data_yield_lookup_output_n2000(ame_table, synthetic_spectrum_path):
    """Larger-N variant matching docs/PERFORMANCE_BASELINE.md's scale.
    Not run by default (see pyproject.toml's pytest markers) -- run
    explicitly with `pytest -m slow` when re-measuring performance."""
    model = _make_model(ame_table, synthetic_spectrum_path, 2000)
    model._calculate_moqs()
    model._calculate_srrf(fref=1.93e6)
    model._simulated_data(harmonics=[127.0], mode="Frequency")
    assert len(model.yield_data) == 2000
    assert np.all(np.asarray(model.yield_data) == 1.0)
```

- [ ] **Step 2: Run it — must pass already (characterization, not new behavior)**

Run: `python3 -m pytest tests/test_analysis.py -v`
Expected: `4 passed` (the `slow`-marked test runs too unless you pass `-m "not slow"`; either way it should pass against the current O(N²) implementation — this task is proving equivalence, not introducing new behavior).

- [ ] **Step 3: Apply the dict-index fix to `_simulated_data`**

Replace:
```python
        self.simulated_data_dict = {}
        self.yield_data = []
        moq_keys = list(self.moq.keys())
        
        for key in moq_keys:
            found = False
            for p in self.particles_to_simulate:
                p_name = f"{int(p[1])}{p[0]}{int(p[4][-1])}+"
                if p_name == key:
                    self.yield_data.append(p[5])
                    found = True
                    break
            if not found: self.yield_data.append(0)
```
with:
```python
        self.simulated_data_dict = {}
        moq_keys = list(self.moq.keys())

        # O(1) name -> yield lookup, built once, instead of an
        # O(len(moq_keys) x len(particles_to_simulate)) nested scan --
        # see docs/PERFORMANCE_BASELINE.md. "First match wins" is
        # preserved deliberately (`if p_name not in yield_by_name`)
        # to match the original loop's `break`-on-first-match semantics
        # exactly, in case particles_to_simulate ever contains a
        # duplicate name.
        yield_by_name = {}
        for p in self.particles_to_simulate:
            p_name = f"{int(p[1])}{p[0]}{int(p[4][-1])}+"
            if p_name not in yield_by_name:
                yield_by_name[p_name] = p[5]
        self.yield_data = [yield_by_name.get(key, 0) for key in moq_keys]
```

- [ ] **Step 4: Re-run — must still pass, unmodified**

Run: `python3 -m pytest tests/test_analysis.py -v`
Expected: `4 passed`, identical assertions, now against the O(N) implementation.

- [ ] **Step 5: Run the full suite**

Run: `python3 -m pytest tests/ -v`
Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add tests/test_analysis.py src/rionid/core.py
git commit -m "Speed fix: O(N^2) -> O(N) yield lookup in _simulated_data

Dict-indexes particles_to_simulate by ion name once instead of scanning
it per moq key. Preserves the original loop's first-match-wins semantics
exactly. tests/test_analysis.py captures identical output at N=10/100
(and N=2000 under -m slow) before and after.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 9: Speed investigation — `plot_all_data` label rendering

**Files:**
- Modify: `src/rionid/gui/plot.py`
- Modify: `tests/test_gui_smoke.py`

**Interfaces:**
- Produces: identical rendered labels (text/position/count), lower per-redraw `QFont` allocation.

This task implements one small, zero-risk fix and explicitly defers a larger one — see rationale below. Do not attempt the larger fix in this task.

- [ ] **Step 1: Extend `tests/test_gui_smoke.py`**

Append:
```python


def test_plot_simulated_data_label_count_and_text(qapp, synthetic_spectrum_path):
    """gui/plot.py's per-label QFont hoist (docs/PERFORMANCE_BASELINE.md)
    must not change any rendered label's text or count."""
    from rionid.core import ImportData
    from rionid.gui.plot import CreatePyGUI
    from rionid.masses import get_ame_data
    from tests.fixtures.synthetic_spectrum import build_ame_candidates

    ame_table = get_ame_data().ame_table
    candidates = build_ame_candidates(ame_table, 10)
    ref_name, ref_aa, ref_zz = candidates[0][0], candidates[0][1], candidates[0][2]
    model = ImportData(
        f"{ref_aa}{ref_name}{ref_zz}+", alphap=0.189,
        filename=synthetic_spectrum_path, reload_data=True,
        circumference=108.36,
    )
    model.ame = get_ame_data()
    model.ame_data = model.ame.ame_table
    model.particles_to_simulate = candidates
    model._calculate_moqs()
    model._calculate_srrf(fref=1.93e6)
    model._simulated_data(harmonics=[127.0], mode="Frequency")

    win = CreatePyGUI()
    win.plot_all_data(model)

    labels = [item.toPlainText() for (_line, item) in win.simulated_items if item is not None]
    assert len(labels) == 10
    assert all(labels)  # every label is non-empty text
```

- [ ] **Step 2: Run it — must pass already**

Run: `python3 -m pytest tests/test_gui_smoke.py -v -k label`
Expected: `1 passed` (characterizes current, correct behavior).

- [ ] **Step 3: Hoist the per-label `QFont` construction in `gui/plot.py`**

In `plot_simulated_data`, this line currently runs inside the per-entry loop, constructing a new `QFont` object for every single label on every redraw:
```python
                text_item = pg.TextItem(text=new_label, color=c, anchor=(0.5, 1))
                text_item.setFont(QFont("Arial", self.font_size))
                text_item.setPos(freq, yield_value * 1.05)
```

Add one line before the `for i, (harmonic, sdata) in enumerate(...)` loop starts:
```python
        color_ref = '#ff7f0e' # Orange for Reference
        color_match = '#2ca02c' # Green for Matches
        label_font = QFont("Arial", self.font_size)
```
(i.e. add `label_font = QFont("Arial", self.font_size)` right after the existing `color_match` line), then inside the per-entry loop replace:
```python
                text_item.setFont(QFont("Arial", self.font_size))
```
with:
```python
                text_item.setFont(label_font)
```

- [ ] **Step 4: Re-run — must still pass**

Run: `python3 -m pytest tests/test_gui_smoke.py -v -k label`
Expected: `1 passed`.

- [ ] **Step 5: Document the deferred, larger optimization**

Add this note to `gui/plot.py`, directly above `plot_simulated_data`, as a code comment (not just in docs — so the next person editing this function sees it immediately):
```python
    def plot_simulated_data(self, data):
        # PERFORMANCE NOTE (see docs/PERFORMANCE_BASELINE.md): profiling
        # at N=2000 candidates shows this method dominated not by
        # pg.TextItem construction itself (~0.14ms/item) but by
        # PyQtGraph's addItem/removeItem scene-graph reparenting overhead
        # (itemChange/changeParent/signal connect-disconnect), triggered
        # on every clear_simulated_data()+plot_simulated_data() cycle --
        # roughly 0.7-1.2s of the ~1.4s total at N=2000. The real fix is
        # reusing TextItem objects across redraws instead of destroying
        # and recreating them, which requires diffing the previous vs.
        # new (harmonic, ion) label set and correctly handling ions that
        # change highlight/reference status between redraws (their
        # colour must update, not just position/text). That correctness
        # cannot be verified visually in this environment (no display
        # server) and was deliberately NOT attempted here -- see
        # docs/superpowers/plans/2026-08-19-wave2a-depid-masses-speed.md
        # Task 9. Revisit once you can visually QA it.
```

- [ ] **Step 6: Run the full suite**

Run: `python3 -m pytest tests/ -v`
Expected: all green.

- [ ] **Step 7: Commit**

```bash
git add src/rionid/gui/plot.py tests/test_gui_smoke.py
git commit -m "Hoist per-label QFont construction out of the redraw loop

Zero-risk, byte-identical-output speed fix: one QFont object is now
reused across all labels in a redraw instead of constructing ~N of them.
Documents (in-code and in the commit) the larger, higher-value TextItem
reuse-across-redraws optimization found by profiling
(itemChange/scene-graph churn dominates at N=2000, ~0.7-1.2s) that was
deliberately NOT attempted -- its correctness under highlight/reference
colour transitions can't be visually verified without a display server
in this environment.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 10: Final re-verification and performance write-up

**Files:**
- Modify: `docs/PERFORMANCE_BASELINE.md`
- Modify: `docs/AUTOMATIC_PID_REMOVAL_MAP.md`

**Interfaces:** None — this task only re-measures and documents; no source changes.

- [ ] **Step 1: Run the full public suite**

Run: `python3 -m pytest tests/ -v`
Expected: all tests pass (`test_correction.py`, `test_masses.py`, `test_gui_smoke.py`, `test_analysis.py`).

Run: `python3 -m pytest tests/ -v -m slow`
Expected: the N=2000 test also passes.

- [ ] **Step 2: Re-run the Phase 0 performance methodology against the now-optimized code**

Adapt and run the same headless timing methodology from `docs/PERFORMANCE_BASELINE.md` (synthetic spectrum, `build_ame_candidates` at N=10/100/2000, `time.perf_counter()` medians) against the current code, covering: `_calculate_moqs`, `_simulated_data`, `AMEData`/`get_ame_data()` construction cost on a second call within the same process, and `plot_all_data` under `QT_QPA_PLATFORM=offscreen`. Record the new numbers.

- [ ] **Step 3: Update `docs/PERFORMANCE_BASELINE.md`**

Add a new section, "Post-Wave-2a re-measurement", with a before/after table for each of the four fixes (moq dict-index, yield dict-index, AME-table + LISEreader double-parse caching, QFont hoist), confirming identical output values (cross-referencing `tests/test_analysis.py`/`test_masses.py`) alongside the new timings. Explicitly note the deferred TextItem-reuse optimization is still open, with a pointer to Task 9's in-code comment.

- [ ] **Step 4: Update `docs/AUTOMATIC_PID_REMOVAL_MAP.md`**

For every row in the "Remove" table, add a short "✅ Removed in Wave 2a" note (or update the table to include a Status column) confirming each item's actual removal, so the document reflects completed work rather than only a plan.

- [ ] **Step 5: Confirm zero automatic-PID / baseline-removal surface remains**

```bash
grep -rniE "quick.?pid|compute_matches|overlay_sim_signal|remove_baseline|psd_baseline" src/ README.md docs/index.md
```
Expected: no matches (aside from this plan/spec/audit documents themselves, which are historical records, not shipped surface).

- [ ] **Step 6: Commit**

```bash
git add docs/PERFORMANCE_BASELINE.md docs/AUTOMATIC_PID_REMOVAL_MAP.md
git commit -m "Wave 2a: re-measure performance, mark automatic-PID removal complete

Confirms the full public suite passes, re-measures the three
performance fixes against docs/PERFORMANCE_BASELINE.md's original
methodology with identical output, and marks every
docs/AUTOMATIC_PID_REMOVAL_MAP.md 'Remove' item as done.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Self-Review Notes (for whoever executes this plan)

- **Spec coverage**: Task 1↔dev tooling; Task 2↔fixtures; Task 3↔correction golden test; Task 4↔masses.py extraction + moq speed fix; Task 5↔PID removal; Task 6↔baseline removal (added during design review); Task 7↔dead I/O; Task 8↔yield speed fix; Task 9↔plot speed investigation; Task 10↔final verification. All six spec scope items and all ten acceptance-criteria checkboxes in the design spec are covered.
- **`_check_io_params`**: intentionally left untouched per your standing instruction — flagged again here so it isn't missed. It becomes dead code the moment Task 5 deletes `quick_pid_script`; deciding its fate is yours, not automated by this plan.
- **`ElBiEn` table fidelity**: Task 4 Step 4's `diff` check is load-bearing — do not skip it. This is physics reference data (electron binding energies), not something to eyeball.
- **Physics-lock**: no task in this plan touches `_calculate_srrf`'s body. `tests/test_correction.py` (Task 3) is the tripwire; every later task's Step "run the full suite" includes it.
