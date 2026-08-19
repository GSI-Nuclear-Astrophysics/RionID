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
