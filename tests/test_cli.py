"""Regression tests for the ``python -m rionid`` command-line interface."""

from pathlib import Path

import pytest

from rionid.__main__ import ESR_CIRCUMFERENCE_M, build_parser, run_controller


def test_cli_requires_candidate_file():
    parser = build_parser()

    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["spectrum.npz", "-r", "72Ge32+", "-ap", "0.51", "-f", "1e6"])

    assert exc_info.value.code == 2


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


def test_run_controller_passes_circumference_and_normalized_frequency_mode(monkeypatch):
    calls = {}

    class FakeData:
        experimental_data = None
        yield_data = []
        nuclei_names = []
        simulated_data_dict = {}
        srrf = []
        ref_frequency = 1.0

        def __init__(self, ref_ion, alphap, filename, circumference):
            calls["init"] = (ref_ion, alphap, filename, circumference)

        def _set_particles_to_simulate_from_file(self, filename):
            calls["particles"] = filename

        def _calculate_moqs(self):
            pass

        def _calculate_srrf(self, **kwargs):
            calls["reference"] = kwargs

        def _simulated_data(self, **kwargs):
            calls["simulation"] = kwargs

    monkeypatch.setattr("rionid.__main__.ImportData", FakeData)

    run_controller(
        "spectrum.npz",
        "ions.lpp",
        0.51,
        "72Ge32+",
        [211.0],
        circumference=128.8,
        fref=1.93e6,
        show=False,
    )

    assert calls["init"] == ("72Ge32+", 0.51, "spectrum.npz", 128.8)
    assert calls["particles"] == "ions.lpp"
    assert calls["simulation"]["mode"] == "frequency"


def test_public_candidate_fixture_runs_full_controller(synthetic_spectrum_path):
    candidate_file = Path(__file__).parents[1] / "examples" / "candidates.lpp"

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
