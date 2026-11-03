"""Regression tests for the ``python -m rionid`` command-line interface."""

from pathlib import Path

from rionid.__main__ import ESR_CIRCUMFERENCE_M, build_parser, run_controller


def test_cli_accepts_manual_candidates_without_lise_file():
    args = build_parser().parse_args(
        [
            "spectrum.npz",
            "-r",
            "72Ge32+",
            "-ap",
            "0.51",
            "-f",
            "1e6",
            "--add-ions",
            "72Ge31+,70Se32+",
        ]
    )
    assert args.filep == ""
    assert args.add_ions == "72Ge31+,70Se32+"


def test_cli_defaults_to_esr_circumference():
    args = build_parser().parse_args(
        [
            "spectrum.npz",
            "-psim",
            "ions.lpp",
            "-r",
            "72Ge32+",
            "-ap",
            "0.51",
            "-f",
            "1e6",
        ]
    )

    assert args.circumference == ESR_CIRCUMFERENCE_M


def test_cli_accepts_custom_circumference():
    args = build_parser().parse_args(
        [
            "spectrum.npz",
            "-psim",
            "ions.lpp",
            "-r",
            "72Ge32+",
            "-ap",
            "0.51",
            "-f",
            "1e6",
            "--circumference",
            "128.8",
        ]
    )

    assert args.circumference == 128.8


def test_run_controller_shares_gui_candidate_pipeline(monkeypatch):
    calls = {}
    sentinel = object()

    def run(**kwargs):
        calls.update(kwargs)
        return sentinel

    monkeypatch.setattr("rionid.__main__.import_controller", run)
    result = run_controller(
        "spectrum.npz",
        "ions.lpp",
        0.51,
        "72Ge32+",
        [211.0],
        circumference=128.8,
        fref=1.93e6,
        show=False,
        manual_ions="72Ge31+",
        display_filter="Ge",
        show_uncertainties=True,
    )
    assert result is sentinel
    assert calls["circumference"] == 128.8
    assert calls["mode"] == "Frequency"
    assert calls["manual_ions"] == "72Ge31+"
    assert calls["display_filter"] == "Ge"
    assert calls["show_uncertainties"]
    assert not calls["write_cache"]


def test_public_candidate_fixture_runs_full_controller(synthetic_spectrum_path):
    candidate_file = Path(__file__).parents[1] / "examples" / "data" / "candidates.lpp"

    model = run_controller(
        synthetic_spectrum_path,
        str(candidate_file),
        0.189,
        "72Ge32+",
        [127.0],
        fref=1.93e6,
        show=False,
    )

    assert set(model.nuclei_names) == {"72Ge32+", "70Se32+", "71Ga31+"}
    assert model.simulated_data_dict["127.0"].shape == (3, 3)


def test_manual_only_full_controller(synthetic_spectrum_path):
    model = run_controller(
        synthetic_spectrum_path,
        "",
        0.189,
        "72Ge32+",
        [127.0],
        fref=1.93e6,
        show=False,
        manual_ions="72Ge31+,70Se32+",
        show_uncertainties=True,
        frequency_scale_hz_per_100kev=50,
        uncertainty_harmonic=127,
    )
    assert set(model.nuclei_names) == {"72Ge32+", "72Ge31+", "70Se32+"}
    assert model.frequency_uncertainties["127.0"]["72Ge31+"] > 0
    assert model.candidate_metadata["72Ge31+"]["source"] == "manual"
