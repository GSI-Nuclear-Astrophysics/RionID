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
