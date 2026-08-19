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
