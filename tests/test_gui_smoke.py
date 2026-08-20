"""GUI regression tests: verify meaningful user-visible state and signal
wiring, not screenshots. Runs headlessly via QT_QPA_PLATFORM=offscreen
(tests/conftest.py sets this before any PyQt5 import).
"""


def test_quick_pid_surface_is_absent(qapp):
    """No automatic-PID controls, methods, or signals may exist."""
    from rionid.gui.inputs import RionID_GUI

    panel = RionID_GUI()

    for attr in (
        "setup_quick_pid",
        "quick_pid_script",
        "onPlotClicked",
        "_stop_quick_pid",
        "overlay_sim_signal",
        "alphap_min_edit",
        "alphap_max_edit",
        "alphap_step_edit",
        "fref_min_edit",
        "fref_max_edit",
        "threshold_edit",
        "matched_result_edit",
    ):
        assert not hasattr(panel, attr), f"automatic-PID attribute still present: {attr}"


def test_collapsible_group_box_is_removed():
    import rionid.gui.inputs as inputs_module

    assert not hasattr(inputs_module, "CollapsibleGroupBox")


def test_compute_matches_is_removed():
    from rionid.core import ImportData

    assert not hasattr(ImportData, "compute_matches")
    assert not hasattr(ImportData, "save_matched_result")


def test_overlay_quick_pid_wiring_is_removed(qapp, qtbot):
    from rionid.gui.app import MainWindow

    win = MainWindow()
    qtbot.addWidget(win)
    assert not hasattr(win, "overlay_simulation")
    assert not hasattr(win.rion_input, "overlay_sim_signal")


def test_visualization_signal_still_updates_the_plot(qapp, qtbot):
    """The retained (non-automatic) run path's signal wiring must survive
    PID removal: RionID_GUI.visualization_signal ->
    MainWindow.update_visualization -> CreatePyGUI.updateData."""
    from rionid.gui.app import MainWindow

    win = MainWindow()
    qtbot.addWidget(win)

    calls = []
    win.visualization_widget.updateData = lambda data: calls.append(data)

    sentinel = object()
    win.rion_input.visualization_signal.emit(sentinel)

    assert calls == [sentinel]


def test_baseline_removal_surface_is_absent(qapp):
    from rionid.gui.inputs import RionID_GUI

    panel = RionID_GUI()
    for attr in ("remove_baseline_checkbox", "psd_baseline_removed_l_edit"):
        assert not hasattr(panel, attr), f"baseline-removal attribute still present: {attr}"


def test_baseline_module_is_removed():
    import importlib

    import pytest

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("rionid.baseline")


def test_peak_detection_surface_is_absent(qapp):
    """Peak detection, its frequency-window controls, and plot markers are removed."""
    from rionid.core import ImportData
    from rionid.gui.inputs import RionID_GUI
    from rionid.gui.plot import CreatePyGUI

    panel = RionID_GUI()
    plot = CreatePyGUI()

    for attr in (
        "peak_thresh_edit",
        "min_distance_edit",
        "matching_freq_min_edit",
        "matching_freq_max_edit",
        "pick_matching_freq_min_button",
        "pick_matching_freq_max_button",
        "enterPlotPickMode",
        "_onPlotPicked",
    ):
        assert not hasattr(panel, attr), f"peak-detection attribute still present: {attr}"

    assert not hasattr(ImportData, "detect_peaks_and_widths")
    assert not hasattr(plot, "plotClicked")
    assert not hasattr(plot, "red_triangles")


def test_font_size_control_is_absent_and_cursor_has_millihertz_precision(qapp):
    from PyQt5.QtWidgets import QLabel, QSpinBox

    from rionid.gui.plot import CreatePyGUI

    plot = CreatePyGUI()

    assert plot.findChildren(QSpinBox) == []
    assert all(label.text() != "Font Size:" for label in plot.findChildren(QLabel))
    assert plot._format_cursor_position(1.930000001) == "Cursor: 1.930000001 MHz"


def test_input_labels_and_fixed_gui_circumference(qapp):
    from PyQt5.QtWidgets import QLabel

    from rionid.core import ImportData
    from rionid.gui.inputs import GUI_RING_CIRCUMFERENCE_M, RionID_GUI

    panel = RionID_GUI()
    labels = {label.text() for label in panel.findChildren(QLabel)}

    assert {
        "Experimental data (.npz):",
        "LISE++ file (.lpp):",
        "Momentum compaction factor (~0.51):",
        "Harmonics (211 212 213):",
        "Reference ion (72Ge32+):",
        "Ions to highlight (72Ge31+ 72Ge30+):",
    } <= labels
    assert not hasattr(panel, "circumference_edit")
    assert GUI_RING_CIRCUMFERENCE_M == 108.36

    model = ImportData("72Ge32+", alphap=0.51, highlight_ions="72Ge31+ 72Ge30+")
    assert model.highlight_ions == ["72Ge31+", "72Ge30+"]


def test_frequency_markers_copy_hz_and_can_be_cleared(qapp):
    from rionid.gui.plot import CreatePyGUI

    plot = CreatePyGUI()
    first = plot.add_frequency_marker(1.930000001)
    second = plot.add_frequency_marker(1.930000002)

    assert first.copy_frequency_hz() == "1930000.001 Hz"
    assert qapp.clipboard().text() == "1930000.001 Hz"
    assert first.frequency_label.toPlainText() == "M1: 1.930000001 MHz"
    assert second.frequency_label.toPlainText() == "M2: 1.930000002 MHz"

    first.setValue(1.930000003)
    assert first.frequency_label.toPlainText() == "M1: 1.930000003 MHz"

    plot.toggle_marker_selection(first)
    plot.toggle_marker_selection(second)
    assert plot.marker_info_label.text() == "Δf M1 ↔ M2: 0.001 Hz"
    assert plot.copy_selected_markers() == "1930000.003 Hz\n1930000.002 Hz"
    assert qapp.clipboard().text() == "1930000.003 Hz\n1930000.002 Hz"

    plot.clear_frequency_markers()
    assert plot.frequency_markers == []
    assert plot.selected_markers == []


def test_action_button_colors(qapp):
    from rionid.gui.inputs import RionID_GUI
    from rionid.gui.plot import CreatePyGUI

    panel = RionID_GUI()
    plot = CreatePyGUI()

    assert "#66bb6a" in panel.run_button.styleSheet()
    assert "#ef5350" in panel.exit_button.styleSheet()
    assert "#bbdefb" in plot.reset_view_button.styleSheet()
    assert "#fff9c4" in plot.clear_markers_button.styleSheet()


def test_plot_simulated_data_label_count_and_text(qapp, synthetic_spectrum_path):
    """The per-label QFont optimization must preserve label text and count."""
    from rionid.core import ImportData
    from rionid.gui.plot import CreatePyGUI
    from rionid.masses import get_ame_data
    from tests.fixtures.synthetic_spectrum import build_ame_candidates

    ame_table = get_ame_data().ame_table
    candidates = build_ame_candidates(ame_table, 10)
    ref_name, ref_aa, ref_zz = candidates[0][0], candidates[0][1], candidates[0][2]
    model = ImportData(
        f"{ref_aa}{ref_name}{ref_zz}+",
        alphap=0.189,
        filename=synthetic_spectrum_path,
        reload_data=True,
        circumference=108.36,
    )
    model.ame = get_ame_data()
    model.ame_data = model.ame.ame_table
    model.particles_to_simulate = candidates
    model._calculate_moqs()
    model._calculate_srrf(fref=1.93e6)
    model._simulated_data(harmonics=[127.0], mode="Frequency")

    win = CreatePyGUI()
    win.plot_all_data(model)

    labels = [item.toPlainText() for (_line, item) in win.simulated_items if item is not None]
    assert len(labels) == 10
    assert all(labels)  # every label is non-empty text
