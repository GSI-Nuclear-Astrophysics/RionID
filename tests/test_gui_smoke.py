"""GUI regression tests: verify meaningful user-visible state and signal
wiring, not screenshots. Runs headlessly via QT_QPA_PLATFORM=offscreen
(tests/conftest.py sets this before any PyQt5 import).
"""


def test_quick_pid_surface_is_absent(qapp):
    """No automatic-PID controls, methods, or signals may exist on the
    input panel -- see docs/AUTOMATIC_PID_REMOVAL_MAP.md."""
    from rionid.gui.inputs import RionID_GUI

    panel = RionID_GUI()

    for attr in (
        "setup_quick_pid", "quick_pid_script", "onPlotClicked",
        "_stop_quick_pid", "overlay_sim_signal",
        "alphap_min_edit", "alphap_max_edit", "alphap_step_edit",
        "fref_min_edit", "fref_max_edit", "threshold_edit",
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


def test_plot_simulated_data_label_count_and_text(qapp, synthetic_spectrum_path):
    """gui/plot.py's per-label QFont hoist (docs/PERFORMANCE_BASELINE.md)
    must not change any rendered label's text or count."""
    from rionid.core import ImportData
    from rionid.gui.plot import CreatePyGUI
    from rionid.masses import get_ame_data
    from tests.fixtures.synthetic_spectrum import build_ame_candidates

    ame_table = get_ame_data().ame_table
    candidates = build_ame_candidates(ame_table, 10)
    ref_name, ref_aa, ref_zz = candidates[0][0], candidates[0][1], candidates[0][2]
    model = ImportData(
        f"{ref_aa}{ref_name}{ref_zz}+", alphap=0.189,
        filename=synthetic_spectrum_path, reload_data=True,
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
