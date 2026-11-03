"""Saved-run regressions using only the redistributable candidate list and invented spectra."""

import csv
import json
from dataclasses import replace
from pathlib import Path

import ezodf
import numpy as np
import pytest

from rionid.bundles import create_run_directory
from rionid.core import ImportData
from rionid.session import (
    RunConfig,
    configuration_values,
    load_configuration,
    replay_configuration,
    resolve_inputs,
    run_configuration,
    save_bundle,
    save_configuration,
)


@pytest.fixture
def run_config(synthetic_spectrum_path, tmp_path):
    # Copy only the public fixture; replay paths never reference unpublished files.
    candidate = tmp_path / "candidates.lpp"
    candidate.write_bytes((Path(__file__).parents[1] / "examples/data/candidates.lpp").read_bytes())
    return RunConfig(
        datafile=synthetic_spectrum_path,
        filep=str(candidate),
        alphap=1.3957,
        refion="72Ge32+",
        harmonics=(125.0, 126.0),
        mode="Frequency",
        value=1.93e6,
        reload_data=True,
        highlight_ions="71Ga31+",
        sim_scalingfactor=2.0,
        nions=2,
        correct=(1e-9, -0.001, 1.0),
        table_filename="my_candidates.out",
    )


def test_completed_bundle_freezes_ionic_inputs_and_replays_without_ame(
    run_config, tmp_path, monkeypatch
):
    completed = run_configuration(run_config)
    folder, files = save_bundle(tmp_path, "public example", run_config, completed)
    assert folder.name.endswith("_RionID_CONFIG_public_example")
    assert {p.name for p in files} == {
        "config.json",
        "candidates.csv",
        "candidates.ods",
        "my_candidates.out",
    }
    payload = json.loads((folder / "config.json").read_text())
    assert set(payload["ionic_inputs"]) == {"72Ge32+", "70Se32+", "71Ga31+"}
    assert payload["configuration"]["circumference"] == 108.36
    monkeypatch.setattr(
        "rionid.core.get_ame_data", lambda: pytest.fail("Replay must use saved masses")
    )
    report = replay_configuration(folder / "config.json")
    assert report.verified
    assert configuration_values(report.completed.config) == configuration_values(run_config)
    assert report.completed.data.simulated_data_dict["126.0"].shape == (2, 3)

    with (folder / "candidates.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6  # all calculated candidates, not only visible top-N
    for row in rows:
        assert float(row["frequency_Hz"]) == pytest.approx(
            float(row["revolution_frequency_Hz"]) * float(row["harmonic"])
        )
    ods = ezodf.opendoc(str(folder / "candidates.ods"))
    assert ods.sheets[0].nrows() == 7
    with (folder / "my_candidates.out").open(newline="") as stream:
        text_rows = list(csv.DictReader(stream, delimiter="\t"))
    assert text_rows == rows


def test_run_does_not_write_source_cache_or_working_directory_output(
    run_config, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    run_configuration(run_config)
    assert not Path("simulation_result.out").exists()
    assert not Path(run_config.datafile).with_name("synthetic_spectrum_cache.npz").exists()


def test_settings_only_roundtrip_and_replay_warning(run_config, tmp_path):
    folder, files = save_bundle(tmp_path, "settings", run_config)
    assert [f.name for f in files] == ["config.json"]
    loaded, payload = load_configuration(files[0])
    assert configuration_values(loaded) == configuration_values(run_config)
    assert "verification" not in payload
    report = replay_configuration(files[0])
    assert not report.verified
    assert report.mismatches == []
    assert any("Settings-only" in w for w in report.warnings)


@pytest.mark.parametrize("which", ["spectrum", "candidates"])
def test_changed_input_rejected_before_calculation(run_config, tmp_path, monkeypatch, which):
    folder, _ = save_bundle(tmp_path, "hash", run_config, run_configuration(run_config))
    payload = json.loads((folder / "config.json").read_text())
    Path(payload["inputs"][which]["path"]).write_bytes(b"changed input")
    monkeypatch.setattr("rionid.session.import_controller", lambda **k: pytest.fail("Must not run"))
    with pytest.raises(ValueError, match="Input hash mismatch"):
        replay_configuration(folder / "config.json")


def test_changed_input_rejected_when_saving_existing_results(run_config, tmp_path):
    completed = run_configuration(run_config)
    Path(run_config.datafile).write_bytes(b"changed")
    with pytest.raises(ValueError, match="Input hash mismatch"):
        save_bundle(tmp_path / "output", "stale", run_config, completed)
    assert not (tmp_path / "output").exists()


def test_cached_spectrum_is_recorded_and_replayed_instead_of_raw(run_config, tmp_path):
    source = Path(run_config.datafile)
    cache = source.with_name(source.stem + "_cache.npz")
    frequencies = np.linspace(9.0, 10.0, 20)
    amplitudes = np.linspace(1.0, 2.0, 20)
    np.savez(cache, frequency=frequencies, amplitude_avg=amplitudes)
    config = replace(run_config, reload_data=False)
    completed = run_configuration(config)
    assert completed.inputs["spectrum"]["path"] == str(cache)
    np.testing.assert_array_equal(completed.data.experimental_data[0], frequencies)
    folder, _ = save_bundle(tmp_path, "cache", config, completed)
    source.unlink()  # the unused source spectrum is not required
    report = replay_configuration(folder / "config.json")
    assert report.verified
    np.testing.assert_array_equal(report.completed.data.experimental_data[0], frequencies)
    assert report.completed.inputs["spectrum"]["sha256"] == completed.inputs["spectrum"]["sha256"]


def test_npz_keys_are_preserved(run_config, tmp_path):
    path = tmp_path / "mapped.npz"
    x, y = np.arange(20.0), np.arange(20.0) + 1
    np.savez(path, custom_frequency=x, custom_amplitude=y)
    config = replace(
        run_config,
        datafile=str(path),
        io_params={"frequency_key": "custom_frequency", "amplitude_key": "custom_amplitude"},
    )
    completed = run_configuration(config)
    folder, _ = save_bundle(tmp_path, "mapped", config, completed)
    report = replay_configuration(folder / "config.json")
    assert report.verified
    np.testing.assert_array_equal(report.completed.data.experimental_data[0], x)


def test_relocated_inputs_match_by_content(run_config, tmp_path):
    folder, _ = save_bundle(tmp_path, "moved", run_config, run_configuration(run_config))
    new_data, new_candidates = tmp_path / "moved.npz", tmp_path / "moved.lpp"
    new_data.write_bytes(Path(run_config.datafile).read_bytes())
    new_candidates.write_bytes(Path(run_config.filep).read_bytes())
    Path(run_config.datafile).unlink()
    Path(run_config.filep).unlink()
    report = replay_configuration(
        folder / "config.json",
        input_overrides={
            "spectrum": str(new_data),
            "candidates": str(new_candidates),
        },
    )
    assert report.verified
    assert report.completed.config.datafile == str(new_data)


def test_modified_checkpoint_is_reported(run_config, tmp_path):
    folder, _ = save_bundle(tmp_path, "check", run_config, run_configuration(run_config))
    config_file = folder / "config.json"
    payload = json.loads(config_file.read_text())
    payload["verification"]["candidates"][0]["frequency_Hz"] += 0.01
    config_file.write_text(json.dumps(payload))
    report = replay_configuration(config_file)
    assert not report.verified
    assert len(report.mismatches) == 1
    assert "frequency_Hz" in report.mismatches[0]


def test_load_does_not_read_inputs(run_config, tmp_path):
    path = tmp_path / "config.json"
    save_configuration(path, run_config)
    Path(run_config.datafile).unlink()
    Path(run_config.filep).unlink()
    loaded, _ = load_configuration(path)
    assert loaded.refion == run_config.refion


def test_bad_schema_and_result_pairing_rejected(run_config, tmp_path):
    path = tmp_path / "config.json"
    save_configuration(path, run_config)
    payload = json.loads(path.read_text())
    payload["software"] = "OtherSoftware"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Unsupported"):
        load_configuration(path)
    completed = run_configuration(run_config)
    with pytest.raises(ValueError, match="different settings"):
        save_configuration(path, replace(run_config, alphap=0.5), completed)


@pytest.mark.parametrize("mode,value", [("Brho", 5.0), ("Kinetic Energy", 350.0), ("Gamma", 1.4)])
def test_all_reference_modes_roundtrip(run_config, tmp_path, mode, value):
    config = replace(run_config, mode=mode, value=value)
    completed = run_configuration(config)
    folder, _ = save_bundle(tmp_path, mode, config, completed)
    report = replay_configuration(folder / "config.json")
    assert report.verified
    assert np.isfinite(completed.data.ref_frequency)


@pytest.mark.parametrize(
    "field,value",
    [
        ("alphap", 0.0),
        ("value", float("nan")),
        ("harmonics", (0.0,)),
        ("nions", -1),
        ("correct", (1.0, 2.0)),
        ("table_filename", "../outside.out"),
        ("table_filename", "config.json"),
    ],
)
def test_invalid_settings_rejected(run_config, field, value):
    with pytest.raises(ValueError):
        replace(run_config, **{field: value}).validate()


def test_unique_safe_folder_names(tmp_path):
    first = create_run_directory(tmp_path, "RionID", "../setting A/")
    second = create_run_directory(tmp_path, "RionID", "../setting A/")
    assert first.parent == tmp_path
    assert first.name.endswith("_RionID_CONFIG_setting_A")
    assert second.name == first.name + "_02"


def test_snapshot_matches_existing_frequency_calculation(run_config):
    # Validate that the new workflow preserves the established scientific calculation.
    old = ImportData(
        run_config.refion,
        1 / run_config.alphap**2,
        filename=run_config.datafile,
        reload_data=True,
        circumference=run_config.circumference,
        write_cache=False,
    )
    old._set_particles_to_simulate_from_file(run_config.filep)
    old._calculate_moqs()
    old._calculate_srrf(fref=run_config.value, correct=list(run_config.correct))
    old._simulated_data(
        harmonics=list(run_config.harmonics),
        mode="Frequency",
        sim_scalingfactor=run_config.sim_scalingfactor,
    )
    new = run_configuration(run_config).data
    np.testing.assert_array_equal(old.srrf, new.srrf)
    assert old.moq == new.moq
    assert old.ref_frequency == new.ref_frequency
    assert resolve_inputs(run_config)["spectrum"]["io_params"]["frequency_key"] == "arr_0"


def test_completed_run_retains_settings_after_caller_edits_config(run_config, tmp_path):
    completed = run_configuration(run_config)
    run_config.alphap = 0.25
    assert completed.config.alphap == 1.3957
    with pytest.raises(ValueError, match="different settings"):
        save_bundle(tmp_path / "wrong", "mismatched", run_config, completed)
    assert not (tmp_path / "wrong").exists()
