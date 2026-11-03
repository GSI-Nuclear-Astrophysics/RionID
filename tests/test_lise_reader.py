"""LISE++ diagnostics preserve overlays while keeping frozen replay offline."""

import csv
import json
import logging

import pytest

from rionid.external.lisereader.reader import LISEreader, MissingMassError
from rionid.session import (
    RunConfig,
    candidate_rows,
    replay_configuration,
    run_configuration,
    save_bundle,
)


def lise_file(tmp_path, records, header="[Calculations]"):
    path = tmp_path / "candidates.lpp"
    path.write_text(header + "\n" + "\n".join(records) + "\n", encoding="latin1")
    return path


def test_duplicate_transport_records_keep_first_match_without_console_noise(
    tmp_path, capsys, caplog
):
    path = lise_file(tmp_path, ["205Au 78 77 =1.0", "205Au 78 76 =9.0"])
    reader = LISEreader(path)
    rows = reader.get_info_all()
    assert rows[0] == rows[1]
    assert rows[0][-2:] == [[78, 77], 1.0]
    assert reader.get_index("205Au78+") == reader.get_index("205Au+78") == 0
    assert reader.summary["usable_records"] == 2
    assert reader.summary["skipped_records"] == 0
    assert capsys.readouterr().out == ""
    assert not caplog.records


def test_missing_masses_emit_one_summary_and_keep_detailed_records(tmp_path, capsys, caplog):
    path = lise_file(
        tmp_path,
        [
            "204Hg 77 =1.0",
            "177Dy 66 =1e-8",
            "178Dy 66 =1e-10",
            "175Tb 65 =1e-10",
            "177Dy 66 =2e-8",
        ],
    )
    reader = LISEreader(path)
    with caplog.at_level(logging.WARNING):
        rows = reader.get_info_all()
    assert len(rows) == 1
    assert rows[0][:4] == ["Hg", 204, 80, 124]
    assert reader.summary == {
        "input_records": 5,
        "usable_records": 1,
        "skipped_records": 4,
        "missing_mass_nuclides": ["175Tb", "177Dy", "178Dy"],
        "invalid_records": 0,
    }
    assert [entry["record"] for entry in reader.skipped_entries] == [2, 3, 4, 5]
    assert {entry["kind"] for entry in reader.skipped_entries} == {"missing_mass"}
    assert len(caplog.records) == 1
    assert "omitted 4 records" in caplog.records[0].getMessage()
    assert capsys.readouterr().out == ""
    with pytest.raises(MissingMassError, match="177Dy"):
        reader.get_info("177Dy66+")


def test_malformed_short_charge_and_yield_records_are_diagnosed(tmp_path, caplog):
    path = lise_file(
        tmp_path,
        [
            "204Hg 77 =1.0",
            "205Au 78",
            "broken",
            "205Au 80 =1.0",
            "206Hg 77 =nan",
            "207Tl 78 =-0.5",
            "208Pb 78 =not-a-number",
        ],
    )
    reader = LISEreader(path)
    with caplog.at_level(logging.WARNING):
        assert len(reader.get_info_all()) == 1
    assert reader.summary["invalid_records"] == 6
    assert reader.summary["missing_mass_nuclides"] == []
    assert len(reader.skipped_entries) == 6
    assert len(caplog.records) == 1
    assert all(entry["reason"] for entry in reader.skipped_entries)


def test_reader_replay_does_not_load_or_query_ame(tmp_path, monkeypatch):
    path = lise_file(tmp_path, ["204Hg 77 =1.0", "177Dy 66 =1e-8"])
    monkeypatch.setattr(
        "rionid.external.lisereader.reader.get_ame_data",
        lambda: pytest.fail("Frozen LISE replay must not query the live AME table"),
    )
    reader = LISEreader(path, validate_masses=False)
    rows = reader.get_info_all()
    assert len(rows) == 2
    assert rows[1][:4] == ["Dy", 177, 66, 111]
    assert reader.skipped_entries == []


def test_skipped_diagnostics_roundtrip_and_export_without_live_ame(
    tmp_path, synthetic_spectrum_path, monkeypatch
):
    lise = lise_file(tmp_path, ["204Hg 77 =1.0", "205Au 77 =0.1", "177Dy 66 =1e-8"])
    config = RunConfig(
        datafile=synthetic_spectrum_path,
        filep=str(lise),
        alphap=0.512,
        refion="204Hg77+",
        harmonics=(213, 214),
        mode="Frequency",
        value=1932523,
        reload_data=True,
    )
    completed = run_configuration(config)
    original_rows = candidate_rows(completed.data)
    folder, files = save_bundle(tmp_path, "diagnostics", config, completed)
    assert "skipped_candidates.csv" in {path.name for path in files}
    payload = json.loads((folder / "config.json").read_text())
    assert payload["lise_import"]["summary"] == completed.data.lise_import_summary
    with (folder / "skipped_candidates.csv").open(newline="") as stream:
        skipped = list(csv.DictReader(stream))
    assert len(skipped) == 1
    assert skipped[0]["ion"] == "177Dy+66"
    assert skipped[0]["kind"] == "missing_mass"
    assert "177Dy66+" not in completed.data.moq
    for module in ("rionid.core", "rionid.masses", "rionid.external.lisereader.reader"):
        monkeypatch.setattr(
            module + ".get_ame_data", lambda: pytest.fail("Frozen replay queried live AME")
        )
    replay = replay_configuration(folder / "config.json")
    assert replay.verified, (replay.warnings, replay.mismatches)
    assert candidate_rows(replay.completed.data) == original_rows
    assert replay.completed.data.lise_import_summary == completed.data.lise_import_summary
    assert replay.completed.data.skipped_lise_candidates == completed.data.skipped_lise_candidates
    second, _ = save_bundle(tmp_path, "replayed", replay.completed.config, replay.completed)
    assert (second / "skipped_candidates.csv").read_bytes() == (
        folder / "skipped_candidates.csv"
    ).read_bytes()


def test_missing_sections_and_empty_calculations_have_clear_errors(tmp_path):
    for content in ("no sections\n", "[Calculations]\n", "[Calculations]\n\n# comment\n"):
        path = tmp_path / "bad.lpp"
        path.write_text(content)
        with pytest.raises(ValueError, match="LISE"):
            LISEreader(path)


def test_only_the_last_calculation_section_is_used(tmp_path):
    path = lise_file(
        tmp_path,
        ["204Hg 77 =1.0", "", "# comment"],
        header=("[Calculations]\nnot a valid record\n[Settings]\nignore this\n[Calculations]"),
    )
    reader = LISEreader(path)
    rows = reader.get_info_all()
    assert len(rows) == 1
    assert rows[0][:4] == ["Hg", 204, 80, 124]


def test_missing_mass_does_not_prevent_manual_candidate_workflow(tmp_path, synthetic_spectrum_path):
    lise = lise_file(tmp_path, ["177Dy 66 =1e-8"])
    config = RunConfig(
        datafile=synthetic_spectrum_path,
        filep=str(lise),
        alphap=0.512,
        refion="204Hg77+",
        manual_ions="205Au76+",
        harmonics=(214,),
        mode="Frequency",
        value=1932523,
        reload_data=True,
    )
    completed = run_configuration(config)
    assert set(completed.data.moq) == {"205Au76+", "204Hg77+"}
    assert completed.data.lise_import_summary["skipped_records"] == 1
    assert completed.data.candidate_metadata["205Au76+"]["source"] == "manual"


def test_intermediate_charge_columns_do_not_change_the_selected_final_ion(tmp_path):
    # Some hand-edited LISE lists retain the original beam's intermediate
    # settings. The final charge is the coordinate used by the overlay.
    path = lise_file(tmp_path, ["159Lu 60 77 77 60 =1.0"])
    reader = LISEreader(path)
    rows = reader.get_info_all()
    assert rows[0][:4] == ["Lu", 159, 71, 88]
    assert rows[0][-2][-1] == 60
    assert reader.skipped_entries == []


def test_diagnostic_export_name_is_reserved(synthetic_spectrum_path):
    config = RunConfig(
        datafile=synthetic_spectrum_path,
        filep="",
        manual_ions="205Au76+",
        refion="204Hg77+",
        harmonics=(214,),
        alphap=0.512,
        mode="Frequency",
        value=1932523,
        table_filename="skipped_candidates.csv",
    )
    with pytest.raises(ValueError, match="Text table name"):
        config.validate()
