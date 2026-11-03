"""Candidate controls, cached display changes and bounded plot complexity."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np

from rionid.gui.inputs import RionID_GUI
from rionid.gui.plot import CreatePyGUI
from rionid.session import RunConfig, configuration_values, run_configuration


def manual_config(path):
    return RunConfig(
        datafile=path,
        filep="",
        alphap=0.189,
        refion="72Ge32+",
        harmonics=(126.0, 127.0),
        mode="Frequency",
        value=1.93e6,
        manual_ions="72Ge31+,70Se32+",
        reload_data=True,
        show_uncertainties=True,
        frequency_scale_hz_per_100kev=25,
        uncertainty_harmonic=127,
    )


def test_extended_form_roundtrip(qapp, synthetic_spectrum_path, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    panel = RionID_GUI()
    config = replace(
        manual_config(synthetic_spectrum_path),
        search_ion="72Ge32+",
        moq_tolerance_percent=0.5,
        charge_min=25,
        charge_max=40,
        display_filter="72Ge",
        highlight_ions="Ge",
        show_isomers=True,
        isomer_filter="205Au, 206Hg",
    )
    panel._restore_run_config(config)
    assert panel.show_isomers_checkbox.isChecked()
    assert panel.isomer_filter_edit.isEnabled()
    assert panel.isomer_filter_edit.text() == "205Au, 206Hg"
    assert configuration_values(panel._current_run_config()) == configuration_values(config)
    assert panel.candidate_options_toggle.isChecked()
    panel.close()


def test_manual_candidates_live_filter_and_mass_bands(
    qapp, synthetic_spectrum_path, monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    completed = run_configuration(manual_config(synthetic_spectrum_path))
    plot = CreatePyGUI()
    panel = RionID_GUI(plot_widget=plot)
    panel._restore_run_config(completed.config)
    panel._accept_completed_run(completed)
    plot.updateData(completed.data)
    assert plot.rendered_uncertainty_count > 0
    panel.display_filter_edit.setText("72Ge")
    panel.highlight_ions_edit.setText("Ge")
    panel._apply_display_options()
    assert plot.rendered_candidate_names == {"72Ge32+", "72Ge31+"}
    assert completed.config.display_filter == "72Ge"
    labels = [item for _, item in plot.simulated_items if item is not None]
    assert any("No LISE++ yield" in item.toolTip() for item in labels)
    panel.display_filter_edit.clear()
    panel.show_uncertainties_checkbox.setChecked(False)
    panel._apply_display_options()
    assert "70Se32+" in plot.rendered_candidate_names
    assert plot.rendered_uncertainty_count == 0
    assert panel.worker is None
    plot.close()
    panel.close()


def test_large_overlay_limits_labels_and_updates_after_zoom(qapp):
    count = 5000
    frequency = np.linspace(405e6, 415e6, count)
    rows = np.column_stack((frequency, np.ones(count), [f"{100 + i}Au75+" for i in range(count)]))
    data = SimpleNamespace(
        experimental_data=(frequency, np.ones(count)),
        simulated_data_dict={"214": rows},
        ref_ion="205Au75+",
        highlight_ions=[],
        candidate_metadata={},
        frequency_uncertainties={},
        display_filter="",
        show_uncertainties=False,
    )
    plot = CreatePyGUI()
    plot.updateData(data)
    assert plot.rendered_candidate_count == count
    assert plot.rendered_label_count <= 300
    plot.plot_widget.setXRange(406, 407, padding=0)
    plot._redraw_candidates()
    assert 400 <= plot.rendered_candidate_count <= 600
    assert plot.rendered_label_count <= 300
    plot.close()


def test_invalid_live_selector_keeps_existing_display(
    qapp, synthetic_spectrum_path, monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    completed = run_configuration(manual_config(synthetic_spectrum_path))
    plot = CreatePyGUI()
    panel = RionID_GUI(plot_widget=plot)
    panel._restore_run_config(completed.config)
    panel._accept_completed_run(completed)
    plot.updateData(completed.data)
    before = set(plot.rendered_candidate_names)
    panel.display_filter_edit.setText("not-an-ion")
    panel._apply_display_options()
    assert plot.rendered_candidate_names == before
    assert completed.config.display_filter == ""
    assert panel.run_status_label.text().startswith("Display unchanged:")
    plot.close()
    panel.close()
