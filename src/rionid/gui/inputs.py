import json
import logging as log
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import toml
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from rionid.candidates import validate_selectors
from rionid.session import (
    RunConfig,
    configuration_values,
    load_configuration,
    replay_configuration,
    resolve_inputs,
    run_configuration,
    save_bundle,
    validate_inputs,
    verification_values,
)

from .dialogs import KeySelectionDialog

log.basicConfig(level=log.DEBUG)

GUI_RING_CIRCUMFERENCE_M = 108.36


class RunWorker(QThread):
    completed = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, operation, parent=None):
        super().__init__(parent)
        self.operation = operation

    def run(self):
        try:
            self.completed.emit(self.operation())
        except Exception as exc:
            self.failed.emit(str(exc))


class RionID_GUI(QWidget):
    """
    The main input control panel for RionID.

    This widget handles file selection, parameter configuration, and the execution
    of simulation scripts. It communicates with the visualization widget to
    handle cursor picking and plot updates.
    """

    visualization_signal = pyqtSignal(object)
    signalError = pyqtSignal(str)

    def __init__(self, plot_widget=None):
        super().__init__()
        self.visualization_widget = plot_widget
        self.saved_data = None
        self.current_io_params = {}
        self._completed_run = None
        self._loaded_config = None
        self._bundle_key = None
        self._bundle_dir = None
        self._bundle_files = []
        self.worker = None
        self._busy = False
        self.initUI()
        self.load_parameters()
        self.signalError.connect(self.show_error)

    def show_error(self, msg):
        QMessageBox.critical(self, "Error", msg)

    def initUI(self):
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        self.vbox = QVBoxLayout(scroll_content)
        self.scroll_area.setWidget(scroll_content)

        main_layout = QVBoxLayout()
        main_layout.addWidget(self.scroll_area)
        self.setLayout(main_layout)

        self.setup_saved_runs()
        self.setup_file_selection()
        self.setup_parameters()
        self.setup_candidate_options()
        self.setup_controls()

    def _form_values(self):
        return {
            "datafile": self.datafile_edit.text(),
            "filep": self.filep_edit.text(),
            "alphap": self.alphap_edit.text(),
            "harmonics": self.harmonics_edit.text(),
            "refion": self.refion_edit.text(),
            "highlight_ions": self.highlight_ions_edit.text(),
            "mode": self.mode_combo.currentText(),
            "value": self.value_edit.text(),
            "sim_scalingfactor": self.sim_scalingfactor_edit.text(),
            "correction": self.correction_edit.text(),
            "nions": self.nions_edit.text(),
            "reload_data": self.reload_data_checkbox.isChecked(),
            "simulation_result": self.simulation_result_edit.text(),
            "io_params": dict(self.current_io_params),
            "output_parent": self.output_parent_edit.text(),
            "folder_label": self.bundle_label_edit.text(),
            "manual_ions": self.manual_ions_edit.text(),
            "search_ion": self.search_ion_edit.text(),
            "moq_tolerance_percent": self.moq_tolerance_percent_edit.text(),
            "charge_min": self.charge_min_edit.text(),
            "charge_max": self.charge_max_edit.text(),
            "display_filter": self.display_filter_edit.text(),
            "frequency_scale_hz_per_100kev": self.frequency_scale_hz_per_100kev_edit.text(),
            "uncertainty_harmonic": self.uncertainty_harmonic_edit.text(),
            "show_uncertainties": self.show_uncertainties_checkbox.isChecked(),
            "show_isomers": self.show_isomers_checkbox.isChecked(),
            "isomer_filter": self.isomer_filter_edit.text(),
        }

    def _restore_form_values(self, values):
        for name in (
            "datafile",
            "filep",
            "alphap",
            "harmonics",
            "refion",
            "highlight_ions",
            "value",
            "sim_scalingfactor",
            "correction",
            "nions",
            "simulation_result",
            "manual_ions",
            "search_ion",
            "display_filter",
            "frequency_scale_hz_per_100kev",
            "charge_max",
        ):
            widget = getattr(self, name + "_edit")
            widget.setText(str(values.get(name, "")))
        mode = values.get("mode", "Frequency")
        self.mode_combo.setCurrentText("Bρ" if mode == "Brho" else mode)
        self.reload_data_checkbox.setChecked(values.get("reload_data", True))
        self.current_io_params = dict(values.get("io_params", {}))
        self.output_parent_edit.setText(values.get("output_parent", ""))
        self.bundle_label_edit.setText(values.get("folder_label", ""))
        for name, default in (
            ("moq_tolerance_percent", "1"),
            ("charge_min", "1"),
            ("uncertainty_harmonic", "1"),
        ):
            getattr(self, name + "_edit").setText(str(values.get(name, default)))
        self.show_uncertainties_checkbox.setChecked(values.get("show_uncertainties", False))
        self.show_isomers_checkbox.setChecked(values.get("show_isomers", False))
        self.isomer_filter_edit.setText(str(values.get("isomer_filter", "")))
        if values.get("manual_ions") or values.get("search_ion") or values.get("show_isomers"):
            self.candidate_options_toggle.setChecked(True)

    def load_parameters(self, filepath="parameters_cache.toml"):
        try:
            with open(filepath) as stream:
                self._restore_form_values(toml.load(stream))
        except FileNotFoundError:
            pass

    def save_parameters(self, filepath="parameters_cache.toml"):
        """Retained for callers using the legacy local settings cache."""
        with open(filepath, "w") as stream:
            toml.dump(self._form_values(), stream)

    def setup_saved_runs(self):
        group = QGroupBox("Saved runs: configuration and tables")
        form = QFormLayout(group)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.output_parent_edit = QLineEdit()
        self.output_parent_edit.setPlaceholderText("Choose an output parent folder")
        self.output_browse_button = QPushButton("Browse…")
        self.output_browse_button.clicked.connect(self.choose_output_parent)
        parent_row = QHBoxLayout()
        parent_row.addWidget(self.output_parent_edit)
        parent_row.addWidget(self.output_browse_button)
        form.addRow("Output parent", parent_row)
        self.bundle_label_edit = QLineEdit()
        self.bundle_label_edit.setPlaceholderText("Optional label, e.g. setting or region")
        form.addRow("Folder label", self.bundle_label_edit)
        self.save_config_button = QPushButton("Save config + tables")
        self.load_config_button = QPushButton("Load config…")
        self.replay_config_button = QPushButton("Replay config…")
        self.save_config_button.clicked.connect(self.save_configuration_bundle)
        self.load_config_button.clicked.connect(self.load_saved_configuration)
        self.replay_config_button.clicked.connect(self.replay_saved_configuration)
        form.addRow(self.save_config_button)
        buttons = QHBoxLayout()
        buttons.addWidget(self.load_config_button)
        buttons.addWidget(self.replay_config_button)
        form.addRow(buttons)
        self.bundle_path_label = QLabel("DATE_RionID_CONFIG[_label]")
        self.bundle_path_label.setWordWrap(True)
        form.addRow(self.bundle_path_label)
        self.vbox.addWidget(group)

    def choose_output_parent(self):
        folder = QFileDialog.getExistingDirectory(self, "Choose output parent folder")
        if folder:
            self.output_parent_edit.setText(folder)

    def setup_file_selection(self):
        self.datafile_label = QLabel("Experimental data (.npz):")
        self.datafile_edit = QLineEdit()
        self.datafile_edit.textChanged.connect(self._reset_io_mapping)
        self.datafile_button = QPushButton("Browse")
        self.datafile_button.clicked.connect(self.browse_datafile)

        self.filep_label = QLabel("LISE++ file (.lpp):")
        self.filep_edit = QLineEdit()
        self.filep_edit.setPlaceholderText("Optional with manual ions or AME search")
        self.filep_button = QPushButton("Browse")
        self.filep_button.clicked.connect(self.browse_lppfile)

        hb1 = QHBoxLayout()
        hb1.addWidget(self.datafile_label)
        hb1.addWidget(self.datafile_edit)
        hb1.addWidget(self.datafile_button)
        self.vbox.addLayout(hb1)

        hb2 = QHBoxLayout()
        hb2.addWidget(self.filep_label)
        hb2.addWidget(self.filep_edit)
        hb2.addWidget(self.filep_button)
        self.vbox.addLayout(hb2)

    def setup_parameters(self):
        # Momentum compaction factor
        self.alphap_edit = QLineEdit()
        hb_ap = QHBoxLayout()
        hb_ap.addWidget(QLabel("Momentum compaction:"))
        hb_ap.addWidget(self.alphap_edit)
        self.vbox.addLayout(hb_ap)

        # Standard Params
        self.harmonics_edit = QLineEdit()
        self.refion_edit = QLineEdit()
        self.highlight_ions_edit = QLineEdit()
        self.highlight_ions_edit.setPlaceholderText("205Au, Au, or 205Au78+")
        self.highlight_ions_edit.setToolTip(
            "An isotope or element selects all its charge states and harmonics."
        )
        self.highlight_ions_edit.editingFinished.connect(self._apply_display_options)

        for lbl, widget in [
            ("Harmonics:", self.harmonics_edit),
            ("Reference ion:", self.refion_edit),
            ("Highlight ions:", self.highlight_ions_edit),
        ]:
            h = QHBoxLayout()
            h.addWidget(QLabel(lbl))
            h.addWidget(widget)
            self.vbox.addLayout(h)

        # Mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Frequency", "Bρ", "Kinetic Energy", "Gamma"])
        self.value_edit = QLineEdit()
        h_mode = QHBoxLayout()
        h_mode.addWidget(QLabel("Mode:"))
        h_mode.addWidget(self.mode_combo)
        h_mode.addWidget(self.value_edit)
        self.vbox.addLayout(h_mode)

        # Scaling Factor
        self.sim_scalingfactor_edit = QLineEdit()
        h_sf = QHBoxLayout()
        h_sf.addWidget(QLabel("Scaling Factor:"))
        h_sf.addWidget(self.sim_scalingfactor_edit)
        self.vbox.addLayout(h_sf)

        # Optional Features Group
        self.optional_group = QGroupBox("Optional Features")
        opt_layout = QFormLayout()

        self.nions_edit = QLineEdit()
        self.nions_edit.setPlaceholderText("e.g. 5")
        opt_layout.addRow("N Ions to Display:", self.nions_edit)

        self.correction_edit = QLineEdit()
        self.correction_edit.setPlaceholderText("a0 a1 a2")
        opt_layout.addRow("Correction (A B C):", self.correction_edit)

        self.reload_data_checkbox = QCheckBox("Reload Data Cache")
        opt_layout.addRow(self.reload_data_checkbox)

        self.simulation_result_edit = QLineEdit()
        self.simulation_result_edit.setPlaceholderText("simulation_result.out")
        opt_layout.addRow("Text table name:", self.simulation_result_edit)

        self.optional_group.setLayout(opt_layout)
        self.vbox.addWidget(self.optional_group)

    def setup_candidate_options(self):
        self.candidate_options_toggle = QToolButton()
        self.candidate_options_toggle.setText("Additional ions, AME search and display")
        self.candidate_options_toggle.setCheckable(True)
        self.candidate_options_toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.candidate_options_toggle.setArrowType(Qt.RightArrow)
        self.vbox.addWidget(self.candidate_options_toggle)
        self.candidate_options = QWidget()
        form = QFormLayout(self.candidate_options)
        form.setRowWrapPolicy(QFormLayout.WrapAllRows)
        self.manual_ions_edit = QLineEdit()
        self.manual_ions_edit.setPlaceholderText("205Au76+, 205Au75+, 201Ir76+")
        self.manual_ions_edit.setToolTip("Add ground-state ions independently of the LISE++ list.")
        form.addRow("Add ions", self.manual_ions_edit)
        self.search_ion_edit = QLineEdit()
        self.search_ion_edit.setPlaceholderText("205Au78+")
        form.addRow("Search AME around m/q of", self.search_ion_edit)
        self.moq_tolerance_percent_edit = QLineEdit("1")
        form.addRow("m/q half-width (%)", self.moq_tolerance_percent_edit)
        self.charge_min_edit = QLineEdit("1")
        self.charge_max_edit = QLineEdit()
        self.charge_max_edit.setPlaceholderText("Z for each nuclide")
        charges = QHBoxLayout()
        charges.addWidget(self.charge_min_edit)
        charges.addWidget(QLabel("to"))
        charges.addWidget(self.charge_max_edit)
        form.addRow("Charge range", charges)
        note = QLabel(
            "Manual and AME candidates have no predicted LISE++ yield; "
            "their marker heights are display guides."
        )
        note.setWordWrap(True)
        form.addRow(note)
        self.show_isomers_checkbox = QCheckBox("Show isomers (Nubase2020)")
        self.show_isomers_checkbox.setToolTip(
            "Adds every listed isomer of each nuclide, in that nuclide's charge states. "
            "Takes effect on the next Run."
        )
        form.addRow(self.show_isomers_checkbox)
        self.isomer_filter_edit = QLineEdit()
        self.isomer_filter_edit.setPlaceholderText("Empty: all nuclides; or e.g. 205Au, 206Hg")
        self.isomer_filter_edit.setToolTip(
            "Only add isomers for these ions, isotopes (205Au) or elements (Au)."
        )
        self.isomer_filter_edit.setEnabled(False)
        self.show_isomers_checkbox.toggled.connect(self.isomer_filter_edit.setEnabled)
        form.addRow("Isomers only for", self.isomer_filter_edit)
        self.display_filter_edit = QLineEdit()
        self.display_filter_edit.setPlaceholderText(
            "Empty: all; 205Au: all its charge states; Au: all gold"
        )
        self.display_filter_edit.editingFinished.connect(self._apply_display_options)
        form.addRow("Show only", self.display_filter_edit)
        self.show_uncertainties_checkbox = QCheckBox("Show ±1σ mass-input frequency bands")
        self.show_uncertainties_checkbox.toggled.connect(self._apply_display_options)
        form.addRow(self.show_uncertainties_checkbox)
        self.frequency_scale_hz_per_100kev_edit = QLineEdit()
        self.frequency_scale_hz_per_100kev_edit.setPlaceholderText(
            "Empty: derivative of corrected frequency model"
        )
        self.frequency_scale_hz_per_100kev_edit.setToolTip(
            "A local approximate conversion, valid in your calibrated region. "
            "The uncertainty of this supplied scale is not propagated."
        )
        form.addRow("Optional Hz per 100 keV", self.frequency_scale_hz_per_100kev_edit)
        self.uncertainty_harmonic_edit = QLineEdit("1")
        self.uncertainty_harmonic_edit.setToolTip(
            "Harmonic at which the Hz/100 keV scale was measured. Use 1 for revolution frequency."
        )
        form.addRow("Scale's harmonic / multiplier", self.uncertainty_harmonic_edit)
        note = QLabel(
            "Bands propagate AME mass inputs. They exclude calibration-model and "
            "electronic-binding uncertainties. Change the scale and press Run to recalculate."
        )
        note.setWordWrap(True)
        form.addRow(note)
        self.candidate_options.hide()
        self.candidate_options_toggle.toggled.connect(self._toggle_candidate_options)
        self.vbox.addWidget(self.candidate_options)

    def _toggle_candidate_options(self, expanded):
        self.candidate_options.setVisible(expanded)
        self.candidate_options_toggle.setArrowType(Qt.DownArrow if expanded else Qt.RightArrow)

    def _apply_display_options(self):
        if self._busy or self.saved_data is None:
            return
        options = dict(
            display_filter=self.display_filter_edit.text().strip(),
            highlight_ions=self.highlight_ions_edit.text().strip(),
            show_uncertainties=self.show_uncertainties_checkbox.isChecked(),
        )
        try:
            validate_selectors(options["display_filter"])
            validate_selectors(options["highlight_ions"])
        except ValueError as exc:
            self.run_status_label.setText(f"Display unchanged: {exc}")
            return
        if self.visualization_widget is not None:
            self.visualization_widget.set_display_options(**options)
        self.saved_data.display_filter = options["display_filter"]
        self.saved_data.highlight_ions = options["highlight_ions"]
        self.saved_data.show_uncertainties = options["show_uncertainties"]
        if self._completed_run is not None:
            self._completed_run.config = replace(self._completed_run.config, **options)

    def setup_controls(self):
        self.run_button = QPushButton("Run")
        self.run_button.setStyleSheet(
            "QPushButton { background-color: #66bb6a; color: black; font-weight: bold; }"
        )
        self.run_button.clicked.connect(self.run_script)
        self.vbox.addWidget(self.run_button)

        self.export_button = QPushButton("Save results bundle (config + tables)")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self.export_results_bundle)
        self.vbox.addWidget(self.export_button)
        self.run_status_label = QLabel("Ready")
        self.run_status_label.setWordWrap(True)
        self.vbox.addWidget(self.run_status_label)

        self.exit_button = QPushButton("Exit")
        self.exit_button.setStyleSheet(
            "QPushButton { background-color: #ef5350; color: black; font-weight: bold; }"
        )
        self.exit_button.clicked.connect(self.close_application)
        self.vbox.addWidget(self.exit_button)

    def close_application(self):
        if not self._busy:
            sys.exit()

    def _get_float(self, widget, default=0.0):
        text = widget.text().strip()
        if not text:
            return default
        try:
            return float(text)
        except ValueError:
            return default

    def _reset_io_mapping(self):
        self.current_io_params = {}

    def _current_run_config(self):
        values = self._form_values()
        correction = values["correction"].strip()
        scaling = values["sim_scalingfactor"].strip()
        nions = values["nions"].strip()
        config = RunConfig(
            datafile=values["datafile"].strip(),
            filep=values["filep"].strip(),
            alphap=float(values["alphap"]),
            refion=values["refion"].strip(),
            harmonics=tuple(float(h) for h in values["harmonics"].replace(",", " ").split()),
            mode="Brho" if values["mode"] == "Bρ" else values["mode"],
            value=float(values["value"]),
            circumference=(
                self._loaded_config.circumference
                if self._loaded_config
                else GUI_RING_CIRCUMFERENCE_M
            ),
            highlight_ions=values["highlight_ions"],
            reload_data=values["reload_data"],
            io_params=dict(self.current_io_params),
            sim_scalingfactor=float(scaling) if scaling else None,
            nions=int(nions) if nions else None,
            correct=tuple(float(c) for c in correction.replace(",", " ").split())
            if correction
            else None,
            table_filename=values["simulation_result"].strip() or "simulation_result.out",
            manual_ions=values["manual_ions"].strip(),
            search_ion=values["search_ion"].strip(),
            moq_tolerance_percent=float(values["moq_tolerance_percent"] or "1"),
            charge_min=int(values["charge_min"] or "1"),
            charge_max=int(values["charge_max"]) if values["charge_max"].strip() else None,
            display_filter=values["display_filter"].strip(),
            frequency_scale_hz_per_100kev=(
                float(values["frequency_scale_hz_per_100kev"])
                if values["frequency_scale_hz_per_100kev"].strip()
                else None
            ),
            uncertainty_harmonic=float(values["uncertainty_harmonic"] or "1"),
            show_uncertainties=values["show_uncertainties"],
            show_isomers=values["show_isomers"],
            isomer_filter=values["isomer_filter"].strip(),
        )
        config.validate()
        return config

    def _check_io_params(self, datafile):
        cache = Path(datafile).with_name(Path(datafile).stem + "_cache.npz")
        if not self.reload_data_checkbox.isChecked() and cache.is_file():
            return True
        if Path(datafile).suffix.lower() == ".npz" and not self.current_io_params:
            with np.load(datafile) as data:
                dlg = KeySelectionDialog(self, list(data.keys()))
            if not dlg.exec_():
                return False
            self.current_io_params = dlg.get_params()
        return True

    def _set_busy(self, busy):
        self._busy = busy
        for button in (
            self.run_button,
            self.save_config_button,
            self.load_config_button,
            self.replay_config_button,
            self.output_browse_button,
        ):
            button.setEnabled(not busy)
        self.export_button.setEnabled(not busy and self._completed_run is not None)

    def _start_operation(self, operation, on_success, status):
        if self._busy:
            return
        self._set_busy(True)
        self.run_status_label.setText(status)
        self.worker = RunWorker(operation, self)
        self.worker.completed.connect(on_success)
        self.worker.failed.connect(self._operation_failed)
        self.worker.finished.connect(self._operation_finished)
        self.worker.start()

    def _operation_finished(self):
        self._set_busy(False)

    def _operation_failed(self, message):
        self.run_status_label.setText("Operation failed")
        self.signalError.emit(message)

    def _accept_completed_run(self, completed):
        self._completed_run = completed
        self.saved_data = completed.data
        self.visualization_signal.emit(completed.data)
        self.run_status_label.setText(
            "Simulation complete. Save config + tables to retain this run."
        )
        skipped = getattr(completed.data, "skipped_lise_candidates", [])
        if skipped:
            self.run_status_label.setText(
                f"Simulation complete. {len(skipped)} LISE++ records were omitted. "
                "Save results for skipped_candidates.csv and the reasons."
            )
            self.run_status_label.setToolTip(
                "\n".join(f"{entry['ion']}: {entry['reason']}" for entry in skipped[:20])
            )
        else:
            self.run_status_label.setToolTip("")

    def run_script(self):
        if self._busy or not self.datafile_edit.text().strip():
            return
        try:
            if not self._check_io_params(self.datafile_edit.text().strip()):
                return
            config = self._current_run_config()
            self._start_operation(
                lambda: run_configuration(config), self._accept_completed_run, "Running simulation…"
            )
        except Exception as exc:
            self.signalError.emit(str(exc))

    def _restore_run_config(self, config):
        parent, label = self.output_parent_edit.text(), self.bundle_label_edit.text()
        self._restore_form_values(
            {
                "datafile": config.datafile,
                "filep": config.filep,
                "alphap": repr(config.alphap),
                "harmonics": " ".join(str(h) for h in config.harmonics),
                "refion": config.refion,
                "highlight_ions": config.highlight_ions,
                "mode": config.mode,
                "value": repr(config.value),
                "sim_scalingfactor": ""
                if config.sim_scalingfactor is None
                else repr(config.sim_scalingfactor),
                "correction": ""
                if config.correct is None
                else " ".join(repr(c) for c in config.correct),
                "nions": "" if config.nions is None else str(config.nions),
                "reload_data": config.reload_data,
                "simulation_result": config.table_filename,
                "io_params": config.io_params,
                "output_parent": parent,
                "folder_label": label,
                "manual_ions": config.manual_ions,
                "search_ion": config.search_ion,
                "moq_tolerance_percent": repr(config.moq_tolerance_percent),
                "charge_min": str(config.charge_min),
                "charge_max": "" if config.charge_max is None else str(config.charge_max),
                "display_filter": config.display_filter,
                "frequency_scale_hz_per_100kev": (
                    ""
                    if config.frequency_scale_hz_per_100kev is None
                    else repr(config.frequency_scale_hz_per_100kev)
                ),
                "uncertainty_harmonic": repr(config.uncertainty_harmonic),
                "show_uncertainties": config.show_uncertainties,
                "show_isomers": config.show_isomers,
                "isomer_filter": config.isomer_filter,
            }
        )
        self._loaded_config = config
        self.candidate_options_toggle.setChecked(
            bool(
                config.manual_ions
                or config.search_ion
                or config.display_filter
                or config.show_uncertainties
                or config.show_isomers
                or config.frequency_scale_hz_per_100kev
            )
        )

    def load_saved_configuration(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load saved configuration", "", "JSON (*.json)")
        if not path:
            return
        try:
            config, _ = load_configuration(path)
            self._restore_run_config(config)
            self.run_status_label.setText(
                "Configuration loaded; press Run or Replay config to calculate."
            )
        except Exception as exc:
            self.signalError.emit(str(exc))

    def replay_saved_configuration(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Replay saved configuration", "", "JSON (*.json)"
        )
        if not path:
            return
        try:
            _, payload = load_configuration(path)
            overrides = {}
            for name, entry in payload["inputs"].items():
                if not Path(entry["path"]).is_file():
                    replacement, _ = QFileDialog.getOpenFileName(
                        self, f"Locate original {name} input", "", "All files (*)"
                    )
                    if not replacement:
                        return
                    overrides[name] = replacement
            self._start_operation(
                lambda: replay_configuration(path, input_overrides=overrides),
                self._accept_replay,
                "Replaying saved run…",
            )
        except Exception as exc:
            self.signalError.emit(str(exc))

    def _accept_replay(self, report):
        self._restore_run_config(report.completed.config)
        self._accept_completed_run(report.completed)
        if report.verified:
            self.run_status_label.setText(
                "Replay verified: inputs and saved numerical results match."
            )
        else:
            messages = report.warnings + report.mismatches
            self.run_status_label.setText("Replay completed with verification notes")
            QMessageBox.warning(self, "Replay verification", "\n".join(messages[:20]))

    def _save_bundle(self, config, completed):
        parent = self.output_parent_edit.text().strip()
        if not parent:
            self.choose_output_parent()
            parent = self.output_parent_edit.text().strip()
        if not parent:
            return
        label = self.bundle_label_edit.text().strip()
        inputs = completed.inputs if completed is not None else resolve_inputs(config)
        validate_inputs(inputs)
        state = {
            "parent": str(Path(parent).expanduser().resolve()),
            "label": label,
            "configuration": configuration_values(config),
            "results": verification_values(completed.data) if completed else None,
            "inputs": inputs,
        }
        key = json.dumps(state, sort_keys=True)
        reusable = key == self._bundle_key
        if not (
            reusable
            and self._bundle_dir is not None
            and all(p.is_file() for p in self._bundle_files)
        ):
            out, files = save_bundle(parent, label, config, completed)
            self._bundle_dir, self._bundle_files, self._bundle_key = out, files, key
        self.bundle_path_label.setText(str(self._bundle_dir))
        self.run_status_label.setText(
            f"Saved configuration and available tables to {self._bundle_dir}"
        )

    def save_configuration_bundle(self):
        if self._busy:
            return
        try:
            if not self._check_io_params(self.datafile_edit.text().strip()):
                return
            config = self._current_run_config()
            completed = self._completed_run
            if completed is not None and configuration_values(config) != configuration_values(
                completed.config
            ):
                completed = None
            self._save_bundle(config, completed)
        except Exception as exc:
            self.signalError.emit(str(exc))

    def export_results_bundle(self):
        if self._busy or self._completed_run is None:
            return
        try:
            self._save_bundle(self._completed_run.config, self._completed_run)
        except Exception as exc:
            self.signalError.emit(str(exc))

    def browse_datafile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Data")
        if f:
            self.datafile_edit.setText(f)

    def browse_lppfile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select LPP")
        if f:
            self.filep_edit.setText(f)

    def closeEvent(self, a0):
        if self._busy:
            a0.ignore()
            return
        super().closeEvent(a0)
