"""Manual candidates, complete charge searches and mass-error propagation."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from rionid.candidates import ion_from_ame, matches_selector, parse_ion, search_moq
from rionid.masses import AMEData, get_ame_data
from rionid.session import (
    RunConfig,
    candidate_rows,
    replay_configuration,
    run_configuration,
    save_bundle,
)


@pytest.fixture
def manual_config(synthetic_spectrum_path):
    return RunConfig(
        datafile=synthetic_spectrum_path,
        filep="",
        alphap=0.512,
        refion="204Hg77+",
        harmonics=(213, 214),
        mode="Frequency",
        value=1932523.3644859814,
        reload_data=True,
        manual_ions="205Au76+, 205Au75+",
        show_uncertainties=True,
    )


def test_ion_and_species_selectors():
    assert parse_ion("205au+76").name == "205Au76+"
    for expression in ("205Au", "205Au*", "Au", "Au*", "205Au76+", "197Os75+,205Au"):
        assert matches_selector("205Au76+", expression)
    assert not matches_selector("205Au76+", "204Au 205Au75+")
    for invalid in ("205Au80+", "205Au0+", "205Au", "2Au1+", "205Xx1+", "205Au76+junk"):
        with pytest.raises(ValueError):
            parse_ion(invalid)


def test_search_includes_all_and_only_interval_and_charge_matches():
    # Compare binary-search results to an independent brute-force enumeration
    # on a small table containing several nuclides and multiple charge states.
    ame = get_ame_data()
    rows = [ame.lookup(e, a) for e, a in [("Hg", 204), ("Au", 205), ("Pt", 204), ("W", 184)]]

    class TinyAME:
        ame_table = rows

        def lookup(self, element, mass):
            return next((r for r in self.ame_table if (r[6], r[5]) == (element, mass)), None)

    table = TinyAME()
    center, _, mass = ion_from_ame("205Au78+", table)
    expected = set()
    for row in rows:
        for charge in range(70, min(row[4], 79) + 1):
            _, _, ionic_mass = ion_from_ame(f"{row[5]}{row[6]}{charge}+", table)
            if abs(ionic_mass / charge / (mass / center.q) - 1) <= 0.01:
                expected.add(f"{row[5]}{row[6]}{charge}+")
    actual = search_moq("205Au78+", 1, 70, 79, table)
    assert {entry["ion"] for entry in actual} == expected
    assert [entry["moq_u_per_e"] for entry in actual] == sorted(
        entry["moq_u_per_e"] for entry in actual
    )
    assert next(entry for entry in actual if entry["ion"] == "205Au78+")["mass_extrapolated"]
    assert search_moq("205Au78+", 0, 78, 78, table)[0]["ion"] == "205Au78+"
    with pytest.raises(ValueError, match="AME2020"):
        search_moq("500Au78+", ame=table)


def test_manual_only_merge_unknown_yields_and_frozen_replay(manual_config, tmp_path, monkeypatch):
    config = replace(manual_config, manual_ions="205Au76+ 205Au75+ 205Au76+", nions=1)
    done = run_configuration(config)
    assert set(done.data.moq) == {"205Au76+", "205Au75+", "204Hg77+"}
    assert len(done.data.simulated_data_dict["214"]) == 3
    for row in candidate_rows(done.data):
        assert row["relative_yield"] is None
        assert row["frequency_sigma_Hz"] >= 0
    assert done.data.frequency_uncertainties["214"]["204Hg77+"] == 0
    folder, _ = save_bundle(tmp_path, "manual", config, done)
    monkeypatch.setattr("rionid.core.get_ame_data", lambda: pytest.fail("Replay requested AME"))
    report = replay_configuration(folder / "config.json")
    assert report.verified
    assert report.completed.config.filep == ""
    assert report.completed.data.candidate_metadata == done.data.candidate_metadata
    assert report.completed.data.frequency_uncertainties == done.data.frequency_uncertainties


def test_lise_duplicate_preserves_yield(manual_config):
    config = replace(
        manual_config,
        refion="72Ge32+",
        manual_ions="72Ge32+ 72Ge31+",
        filep=str(Path(__file__).parents[1] / "examples/data/candidates.lpp"),
    )
    done = run_configuration(config)
    metadata = done.data.candidate_metadata["72Ge32+"]
    assert metadata["source"] == "lise"
    assert metadata["sources"] == ["lise", "manual"]
    assert metadata["relative_yield_known"]
    assert not done.data.candidate_metadata["72Ge31+"]["relative_yield_known"]


def test_reference_is_auto_added_outside_search_interval(manual_config):
    config = replace(
        manual_config, manual_ions="", refion="1H1+", search_ion="205Au78+", charge_min=75
    )
    done = run_configuration(config)
    assert done.data.candidate_metadata["1H1+"]["source"] == "reference"
    assert all(abs(entry["offset_percent"]) <= 1 for entry in done.data.search_results)


def test_uncertainty_includes_calibration_derivative_and_reference_mass(manual_config):
    cfg = replace(manual_config, correct=(2.6794364554e-7, -1.03579554, 1001024.3598130841))
    model = run_configuration(cfg).data
    target = "205Au76+"
    mu, ref = model.moq[target], model.moq[model.ref_ion]
    target_sigma = model.candidate_metadata[target]["mass_uncertainty_keV"] / (
        76 * AMEData.UU * 1000
    )
    ref_sigma = model.candidate_metadata[model.ref_ion]["mass_uncertainty_keV"] / (
        77 * AMEData.UU * 1000
    )
    f0 = cfg.value * (1 - cfg.alphap * (mu / ref - 1))
    slope = 1 + 2 * cfg.correct[0] * f0 + cfg.correct[1]
    expected = abs(slope) * np.hypot(
        cfg.alphap * cfg.value / ref * target_sigma,
        cfg.alphap * cfg.value * mu / ref**2 * ref_sigma,
    )
    assert model.frequency_uncertainties["214"][target] == pytest.approx(214 * expected)


def test_charge_states_share_their_atomic_mass_error(manual_config):
    cfg = replace(manual_config, refion="205Au77+", manual_ions="205Au76+")
    data = run_configuration(cfg).data
    mu, reference = data.moq["205Au76+"], data.moq["205Au77+"]
    sigma = data.candidate_metadata["205Au76+"]["mass_uncertainty_keV"]
    target_term = -cfg.alphap * cfg.value / reference * sigma / (76 * AMEData.UU * 1000)
    ref_term = cfg.alphap * cfg.value * mu / reference**2 * sigma / (77 * AMEData.UU * 1000)
    expected = 214 * abs(target_term + ref_term)
    assert data.frequency_uncertainties["214"]["205Au76+"] == pytest.approx(expected)
    assert expected < 0.01


def test_empirical_uncertainty_scale_has_explicit_harmonic(manual_config):
    cfg = replace(manual_config, frequency_scale_hz_per_100kev=10, uncertainty_harmonic=214)
    model = run_configuration(cfg).data
    assert model.frequency_uncertainties["214"]["205Au76+"] == pytest.approx(20)
    assert model.frequency_uncertainties["213"]["205Au76+"] == pytest.approx(20 * 213 / 214)


def test_rf_and_revolution_coordinates_give_identical_uncertainty(manual_config):
    fundamental = replace(manual_config, correct=(2.6794364554e-7, -1.03579554, 1001024.3598130841))
    rf = replace(
        manual_config,
        value=413560000,
        harmonics=(213 / 214, 1),
        correct=(1.25207311e-9, -1.03579554, 2.14219213e8),
    )
    first, second = run_configuration(fundamental).data, run_configuration(rf).data
    for h1, h2 in zip(fundamental.harmonics, rf.harmonics):
        for ion in first.moq:
            assert first.frequency_uncertainties[str(h1)][ion] == pytest.approx(
                second.frequency_uncertainties[str(h2)][ion], abs=1e-9
            )


def test_rigidity_does_not_depend_on_overlay_harmonic(manual_config):
    first = run_configuration(replace(manual_config, harmonics=(1,))).data
    second = run_configuration(replace(manual_config, harmonics=(214,))).data
    assert first.brho == pytest.approx(second.brho)
    assert first.brho > 5


@pytest.mark.parametrize(
    "changes",
    [
        {"manual_ions": "205Au80+"},
        {"search_ion": "205Au"},
        {"moq_tolerance_percent": -1},
        {"charge_min": 0},
        {"charge_max": 0},
        {"frequency_scale_hz_per_100kev": -1},
        {"uncertainty_harmonic": 0},
        {"display_filter": "invalidstring"},
    ],
)
def test_candidate_config_validation(manual_config, changes):
    with pytest.raises(ValueError):
        replace(manual_config, **changes).validate()


@pytest.mark.parametrize(
    "mode,value",
    [
        ("Frequency", 1932523.364),
        ("Brho", 7.8),
        ("Kinetic Energy", 400),
        ("Gamma", 1.4),
    ],
)
def test_reference_mode_derivatives_against_mass_perturbations(manual_config, mode, value):
    from copy import deepcopy

    cfg = replace(manual_config, mode=mode, value=value, alphap=0.3, correct=(2e-7, -0.7, 6e5))
    model = run_configuration(cfg).data
    target = "205Au76+"
    index = list(model.moq).index(target)
    kwargs = {
        "Frequency": {"fref": value},
        "Brho": {"brho": value},
        "Kinetic Energy": {"ke": value},
        "Gamma": {"gam": value},
    }[mode]
    terms = []
    step_keV = 5
    for changed in (target, cfg.refion):
        charge = parse_ion(changed).q
        values = []
        for direction in (-1, 1):
            shifted = deepcopy(model)
            shifted.moq[changed] += direction * step_keV / (charge * AMEData.UU * 1000)
            shifted._calculate_srrf(**kwargs, correct=list(cfg.correct))
            values.append(shifted.srrf[index] * shifted.ref_frequency)
        derivative = (values[1] - values[0]) / (2 * step_keV)
        terms.append(derivative * model.candidate_metadata[changed]["mass_uncertainty_keV"])
    expected = np.hypot(*terms) * 214
    assert model.frequency_uncertainties["214"][target] == pytest.approx(expected, rel=3e-7)


def test_reference_alternate_spelling_is_canonical(manual_config):
    data = run_configuration(replace(manual_config, refion="204hg+77")).data
    assert data.ref_ion == "204Hg77+"
    assert data.frequency_uncertainties["214"][data.ref_ion] == 0
