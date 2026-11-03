"""GUI saving/replay using public candidates and invented spectra."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from rionid.gui.inputs import RionID_GUI
from rionid.session import RunConfig, configuration_values, run_configuration, save_bundle


@pytest.fixture
def run_config(synthetic_spectrum_path):
    return RunConfig(
        datafile=synthetic_spectrum_path,
        filep=str(Path(__file__).parents[1] / "examples/data/candidates.lpp"),
        alphap=0.51,
        refion="72Ge32+",
        harmonics=(125.0, 126.0),
        mode="Frequency",
        value=1.93e6,
        reload_data=True,
        io_params={"frequency_key": "arr_0", "amplitude_key": "arr_1"},
        correct=(1e-9, -0.001, 1.0),
        highlight_ions="71Ga31+",
        nions=2,
    )


@pytest.fixture
def panel(qtbot, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)  # no real local settings cache participates in tests
    panel = RionID_GUI()
    qtbot.addWidget(panel)
    panel.output_parent_edit.setText(str(tmp_path / "saved_runs"))
    panel.bundle_label_edit.setText("region A")
    return panel


def wait_finished(panel, qtbot):
    qtbot.waitUntil(lambda: not panel._busy, timeout=10000)
    assert panel.worker.wait(1000)


def test_visible_controls_in_main_window(qtbot, monkeypatch, tmp_path, qapp):
    from rionid.gui.app import MainWindow

    monkeypatch.chdir(tmp_path)
    win = MainWindow()
    qtbot.addWidget(win)
    win.resize(1400, 900)
    win.show()
    qapp.processEvents()
    viewport = win.rion_input.scroll_area.viewport()
    for button in (
        win.rion_input.save_config_button,
        win.rion_input.replay_config_button,
        win.rion_input.load_config_button,
        win.rion_input.output_browse_button,
    ):
        assert viewport.rect().contains(button.mapTo(viewport, button.rect().topLeft()))
        assert viewport.rect().contains(button.mapTo(viewport, button.rect().bottomRight()))


def test_restore_preserves_all_config_settings_and_output_parent(panel, run_config):
    config = replace(
        run_config,
        circumference=128.8,
        mode="Gamma",
        value=1.4,
        sim_scalingfactor=2.5,
        table_filename="custom.out",
    )
    output_parent = panel.output_parent_edit.text()
    panel._restore_run_config(config)
    assert configuration_values(panel._current_run_config()) == configuration_values(config)
    assert panel.output_parent_edit.text() == output_parent
    assert panel.bundle_label_edit.text() == "region A"


def test_run_save_and_export_buttons_share_dated_folder(panel, run_config, qtbot):
    panel._restore_run_config(run_config)
    qtbot.mouseClick(panel.run_button, Qt.LeftButton)
    wait_finished(panel, qtbot)
    assert panel.saved_data is not None
    assert panel.export_button.isEnabled()
    qtbot.mouseClick(panel.save_config_button, Qt.LeftButton)
    folder = panel._bundle_dir
    assert folder.name.endswith("_RionID_CONFIG_region_A")
    assert {p.name for p in folder.iterdir()} == {
        "config.json",
        "candidates.csv",
        "candidates.ods",
        "simulation_result.out",
    }
    qtbot.mouseClick(panel.export_button, Qt.LeftButton)
    qtbot.mouseClick(panel.save_config_button, Qt.LeftButton)
    assert panel._bundle_dir == folder
    assert len(list(folder.parent.iterdir())) == 1


def test_edited_settings_save_without_old_results_and_export_keeps_original(
    panel, run_config, qtbot
):
    panel._restore_run_config(run_config)
    panel._accept_completed_run(run_configuration(run_config))
    panel.alphap_edit.setText("0.5")
    qtbot.mouseClick(panel.save_config_button, Qt.LeftButton)
    config_only = panel._bundle_dir
    payload = json.loads((config_only / "config.json").read_text())
    assert payload["configuration"]["alphap"] == 0.5
    assert "verification" not in payload
    qtbot.mouseClick(panel.save_config_button, Qt.LeftButton)
    assert panel._bundle_dir == config_only  # repeated settings-only save also reuses its folder
    panel.export_results_bundle()
    results_folder = panel._bundle_dir
    payload = json.loads((results_folder / "config.json").read_text())
    assert results_folder != config_only
    assert payload["configuration"]["alphap"] == 0.51
    assert "verification" in payload
    assert (results_folder / "candidates.csv").exists()


def test_replay_button_runs_restores_and_verifies(panel, run_config, tmp_path, qtbot, monkeypatch):
    folder, _ = save_bundle(tmp_path, "original", run_config, run_configuration(run_config))
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", lambda *a, **k: (str(folder / "config.json"), "")
    )
    displayed = []
    panel.visualization_signal.connect(displayed.append)
    qtbot.mouseClick(panel.replay_config_button, Qt.LeftButton)
    wait_finished(panel, qtbot)
    assert panel.run_status_label.text().startswith("Replay verified")
    assert configuration_values(panel._current_run_config()) == configuration_values(run_config)
    assert len(displayed) == 1
    assert displayed[0] is panel.saved_data


def test_load_config_restores_without_calculation(panel, run_config, tmp_path, qtbot, monkeypatch):
    folder, _ = save_bundle(tmp_path, "settings", run_config)
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", lambda *a, **k: (str(folder / "config.json"), "")
    )
    qtbot.mouseClick(panel.load_config_button, Qt.LeftButton)
    assert panel.worker is None
    assert panel.saved_data is None
    assert configuration_values(panel._current_run_config()) == configuration_values(run_config)


def test_replay_failure_recovers_controls(panel, run_config, tmp_path, qtbot, monkeypatch):
    folder, _ = save_bundle(tmp_path, "original", run_config, run_configuration(run_config))
    Path(run_config.datafile).write_bytes(b"changed spectrum")
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", lambda *a, **k: (str(folder / "config.json"), "")
    )
    errors = []
    monkeypatch.setattr(QMessageBox, "critical", lambda *args: errors.append(args[-1]))
    qtbot.mouseClick(panel.replay_config_button, Qt.LeftButton)
    wait_finished(panel, qtbot)
    assert errors and "Input hash mismatch" in errors[0]
    assert panel.run_button.isEnabled()
    assert panel.save_config_button.isEnabled()
    assert panel.replay_config_button.isEnabled()
    assert panel.saved_data is None


def test_data_path_change_resets_npz_mapping(panel, run_config):
    panel._restore_run_config(run_config)
    assert panel.current_io_params
    panel.datafile_edit.setText("different.npz")
    assert panel.current_io_params == {}


def test_legacy_parameter_cache_roundtrip_preserves_key_mapping(panel, run_config, tmp_path):
    panel._restore_run_config(run_config)
    path = tmp_path / "legacy.toml"
    panel.save_parameters(path)
    panel.alphap_edit.setText("0.4")
    panel.load_parameters(path)
    assert configuration_values(panel._current_run_config()) == configuration_values(run_config)


def test_busy_main_window_cannot_destroy_running_thread(qtbot, monkeypatch, tmp_path):
    from rionid.gui.app import MainWindow

    monkeypatch.chdir(tmp_path)
    win = MainWindow()
    qtbot.addWidget(win)
    win.show()
    win.rion_input._set_busy(True)
    assert not win.close()
    assert win.isVisible()
    win.rion_input._set_busy(False)
    assert win.close()
