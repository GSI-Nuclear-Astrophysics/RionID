"""Nubase2020 isomers: parsing, ion names and the Show isomers option."""

from dataclasses import replace

import pytest

from rionid.__main__ import run_controller
from rionid.candidates import ion_from_ame, matches_selector, parse_ion
from rionid.isomers import NubaseData
from rionid.masses import AMEData
from rionid.session import RunConfig, replay_configuration, run_configuration, save_bundle

NUBASE_LINES = """\
# header line
072 0320   72Ge   -72585.91       0.08                                stbl              0+
072 0321   72Ge m -71894.48       0.09      691.43       0.04         444.2   ns 0.8    0+
180 0721   180Hfm -48637.9        1.4      1141.552      0.015          5.53   h 0.02   8-
180 0722   180Hfn -48405.1        1.4      1374.36       0.04         570     us 20     4-
089 0391   89Y  m -86802.2        0.3       908.97       0.03          15.663  s 0.005  9/2+*
089 0411   89Nb m -80630#        40#          0#        30#        *    1.10   h 0.03   (1/2)-
016 0083   16O  p   8231.60       0.27    12968.6        0.27                           2-
090 0413   90Nb p -82491          3         171.10       0.10          <1     us        7+
"""


@pytest.fixture
def nubase_file(tmp_path):
    path = tmp_path / "nubase.data"
    path.write_text(NUBASE_LINES)
    return str(path)


@pytest.fixture
def ge_isomer_config(synthetic_spectrum_path):
    return RunConfig(
        datafile=synthetic_spectrum_path,
        filep="",
        alphap=0.189,
        refion="72Ge32+",
        harmonics=(127.0,),
        mode="Frequency",
        value=1.93e6,
        reload_data=True,
        manual_ions="72Ge32+ 72Ge31+",
        show_isomers=True,
    )


def test_parse_keeps_isomers_and_skips_levels(nubase_file):
    nubase = NubaseData(nubase_file)
    (ge,) = nubase.isomers("Ge", 72)
    assert (ge.label, ge.excitation_keV, ge.excitation_unc_keV) == ("m", 691.43, 0.04)
    assert [i.label for i in nubase.isomers("Hf", 180)] == ["m", "n"]
    assert nubase.isomers("O", 16) == []
    assert [i.label for i in nubase.isomers("Nb", 90)] == ["p"]
    assert nubase.isomers("Nb", 89) == []  # excitation energy 0# is not a line position
    assert not ge.extrapolated


def test_isomer_ion_names():
    ion = parse_ion("180Hfm70+")
    assert (ion.A, ion.element, ion.isomer, ion.q, ion.name) == (180, "Hf", "m", 70, "180Hfm70+")
    assert parse_ion("96Cm20+").isomer == "" and parse_ion("96Cm20+").element == "Cm"
    assert parse_ion("180hfm+70").name == "180Hfm70+"
    assert matches_selector("180Hfm70+", "180Hf") and matches_selector("180Hfm70+", "Hf")
    assert not matches_selector("180Hfm70+", "180Hf70+")
    for invalid in ("180Xxm76+", "180Hfm80+", "180Hfmm70+"):
        with pytest.raises(ValueError):
            parse_ion(invalid)
    with pytest.raises(ValueError, match="isomer"):
        ion_from_ame("180Hfm70+")


def test_isomers_added_for_each_charge_state(ge_isomer_config, monkeypatch, nubase_file):
    monkeypatch.setattr("rionid.core.get_nubase_data", lambda: NubaseData(nubase_file))
    data = run_configuration(ge_isomer_config).data
    assert set(data.moq) == {"72Ge32+", "72Ge31+", "72Gem32+", "72Gem31+"}
    for q in (31, 32):
        gap = data.total_mass[f"72Gem{q}+"] - data.total_mass[f"72Ge{q}+"]
        assert gap == pytest.approx(AMEData.to_u(0.69143), rel=1e-9)
        assert data.moq[f"72Gem{q}+"] == pytest.approx(data.total_mass[f"72Gem{q}+"] / q)
    metadata = data.candidate_metadata["72Gem32+"]
    assert metadata["source"] == "isomer" and not metadata["relative_yield_known"]
    assert metadata["mass_uncertainty_keV"] > 0.04
    assert data.protons["72Gem32+"] == 32
    assert len(data.simulated_data_dict["127.0"]) == 4


def test_isomer_filter_limits_nuclides_and_charge_states(
    ge_isomer_config, monkeypatch, nubase_file
):
    monkeypatch.setattr("rionid.core.get_nubase_data", lambda: NubaseData(nubase_file))
    config = replace(ge_isomer_config, manual_ions="72Ge32+ 72Ge31+ 180Hf70+")

    everything = set(run_configuration(config).data.moq)
    assert {"72Gem32+", "72Gem31+", "180Hfm70+", "180Hfn70+"} <= everything

    only_hf = set(run_configuration(replace(config, isomer_filter="180Hf, 206Hg")).data.moq)
    assert only_hf == {"72Ge32+", "72Ge31+", "180Hf70+", "180Hfm70+", "180Hfn70+"}

    one_charge = set(run_configuration(replace(config, isomer_filter="Ge")).data.moq)
    assert one_charge - {"72Ge32+", "72Ge31+", "180Hf70+"} == {"72Gem32+", "72Gem31+"}

    ion = set(run_configuration(replace(config, isomer_filter="72Ge31+")).data.moq)
    assert ion - {"72Ge32+", "72Ge31+", "180Hf70+"} == {"72Gem31+"}

    with pytest.raises(ValueError):
        replace(config, isomer_filter="180Hf76+junk").validate()


def test_isomers_off_by_default(ge_isomer_config, monkeypatch):
    monkeypatch.setattr(
        "rionid.core.get_nubase_data", lambda: pytest.fail("Nubase read without the option")
    )
    data = run_configuration(replace(ge_isomer_config, show_isomers=False)).data
    assert set(data.moq) == {"72Ge32+", "72Ge31+"}


def test_isomer_run_replays_without_nubase(ge_isomer_config, tmp_path, monkeypatch, nubase_file):
    monkeypatch.setattr("rionid.core.get_nubase_data", lambda: NubaseData(nubase_file))
    done = run_configuration(ge_isomer_config)
    folder, _ = save_bundle(tmp_path, "isomers", ge_isomer_config, done)
    monkeypatch.setattr(
        "rionid.core.get_nubase_data", lambda: pytest.fail("Replay requested Nubase")
    )
    report = replay_configuration(folder / "config.json")
    assert report.verified
    assert report.completed.config.show_isomers
    assert set(report.completed.data.moq) == set(done.data.moq)


def test_cli_controller_forwards_show_isomers(monkeypatch):
    calls = {}
    monkeypatch.setattr("rionid.__main__.import_controller", lambda **kw: calls.update(kw))
    run_controller(
        "s.npz",
        "",
        0.5,
        "72Ge32+",
        [1.0],
        fref=1.9e6,
        show=False,
        show_isomers=True,
        isomer_filter="205Au",
    )
    assert calls["show_isomers"] and calls["isomer_filter"] == "205Au"
