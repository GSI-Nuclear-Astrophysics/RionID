import argparse
import logging as log
import os
import sys

import numpy as np
import toml
from PyQt5.QtCore import pyqtSignal
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
    QVBoxLayout,
    QWidget,
)

from .controller import import_controller
from .dialogs import KeySelectionDialog

log.basicConfig(level=log.DEBUG)

GUI_RING_CIRCUMFERENCE_M = 108.36


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

        self.setup_file_selection()
        self.setup_parameters()
        self.setup_controls()

    def load_parameters(self, filepath="parameters_cache.toml"):
        try:
            with open(filepath, "r") as f:
                p = toml.load(f)
                self.datafile_edit.setText(p.get("datafile", ""))
                self.filep_edit.setText(p.get("filep", ""))
                self.alphap_edit.setText(p.get("alphap", ""))
                self.harmonics_edit.setText(p.get("harmonics", ""))
                self.refion_edit.setText(p.get("refion", ""))
                self.highlight_ions_edit.setText(p.get("highlight_ions", ""))
                self.mode_combo.setCurrentText(p.get("mode", "Frequency"))
                self.value_edit.setText(p.get("value", ""))
                self.sim_scalingfactor_edit.setText(p.get("sim_scalingfactor", ""))
                self.correction_edit.setText(p.get("correction", ""))
                self.nions_edit.setText(p.get("nions", ""))
                self.reload_data_checkbox.setChecked(p.get("reload_data", True))
                self.simulation_result_edit.setText(p.get("simulation_result", ""))
        except FileNotFoundError:
            pass

    def save_parameters(self, filepath="parameters_cache.toml"):
        p = {
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
        }
        with open(filepath, "w") as f:
            toml.dump(p, f)

    def setup_file_selection(self):
        self.datafile_label = QLabel("Experimental data (.npz):")
        self.datafile_edit = QLineEdit()
        self.datafile_button = QPushButton("Browse")
        self.datafile_button.clicked.connect(self.browse_datafile)

        self.filep_label = QLabel("LISE++ file (.lpp):")
        self.filep_edit = QLineEdit()
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
        hb_ap.addWidget(QLabel("Momentum compaction factor (~0.51):"))
        hb_ap.addWidget(self.alphap_edit)
        self.vbox.addLayout(hb_ap)

        # Standard Params
        self.harmonics_edit = QLineEdit()
        self.refion_edit = QLineEdit()
        self.highlight_ions_edit = QLineEdit()

        for lbl, widget in [
            ("Harmonics (211 212 213):", self.harmonics_edit),
            ("Reference ion (72Ge32+):", self.refion_edit),
            ("Ions to highlight (72Ge31+ 72Ge30+):", self.highlight_ions_edit),
        ]:
            h = QHBoxLayout()
            h.addWidget(QLabel(lbl))
            h.addWidget(widget)
            self.vbox.addLayout(h)

        # Mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Frequency", "Bρ", "Kinetic Energy"])
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
        opt_layout.addRow("Correction (a0*x**2 + a1*x + a2):", self.correction_edit)

        self.reload_data_checkbox = QCheckBox("Reload Data Cache")
        opt_layout.addRow(self.reload_data_checkbox)

        self.simulation_result_edit = QLineEdit()
        opt_layout.addRow("Sim Result File:", self.simulation_result_edit)

        self.optional_group.setLayout(opt_layout)
        self.vbox.addWidget(self.optional_group)

    def setup_controls(self):
        self.run_button = QPushButton("Run")
        self.run_button.setStyleSheet(
            "QPushButton { background-color: #66bb6a; color: black; font-weight: bold; }"
        )
        self.run_button.clicked.connect(self.run_script)
        self.vbox.addWidget(self.run_button)

        self.exit_button = QPushButton("Exit")
        self.exit_button.setStyleSheet(
            "QPushButton { background-color: #ef5350; color: black; font-weight: bold; }"
        )
        self.exit_button.clicked.connect(self.close_application)
        self.vbox.addWidget(self.exit_button)

    def close_application(self):
        sys.exit()

    def _get_float(self, widget, default=0.0):
        text = widget.text().strip()
        if not text:
            return default
        try:
            return float(text)
        except ValueError:
            return default

    def run_script(self):
        datafile = self.datafile_edit.text()
        if not datafile:
            return

        io_params = self.current_io_params
        ext = os.path.splitext(datafile)[1].lower()
        if ext == ".npz" and not io_params:
            data = np.load(datafile)
            dlg = KeySelectionDialog(self, list(data.keys()))
            if dlg.exec_():
                io_params = dlg.get_params()
                self.current_io_params = io_params
            else:
                return

        correct_str = self.correction_edit.text().strip()
        correct = [float(x) for x in correct_str.split()] if correct_str else None

        sim_sf_str = self.sim_scalingfactor_edit.text().strip()
        sim_sf = float(sim_sf_str) if sim_sf_str else None

        alphap = self._get_float(self.alphap_edit, 0.0)

        args = argparse.Namespace(
            datafile=datafile,
            filep=self.filep_edit.text(),
            alphap=alphap,
            harmonics=self.harmonics_edit.text(),
            refion=self.refion_edit.text(),
            circumference=GUI_RING_CIRCUMFERENCE_M,
            mode=self.mode_combo.currentText(),
            value=self.value_edit.text(),
            highlight_ions=self.highlight_ions_edit.text(),
            io_params=io_params,
            reload_data=self.reload_data_checkbox.isChecked(),
            nions=self.nions_edit.text(),
            sim_scalingfactor=sim_sf,
            correct=correct,
        )

        try:
            self.save_parameters()
            data = import_controller(**vars(args))
            self.saved_data = data
            self.visualization_signal.emit(data)
        except Exception as e:
            self.signalError.emit(str(e))

    def _check_io_params(self, datafile):
        """
        Helper method to ensure NPZ files have their keys mapped before running scripts.
        Returns True if parameters are set or not needed, False if the user cancels.
        """
        ext = os.path.splitext(datafile)[1].lower()
        if ext == ".npz" and not self.current_io_params:
            try:
                data = np.load(datafile)
                dlg = KeySelectionDialog(self, list(data.keys()))
                if dlg.exec_():
                    self.current_io_params = dlg.get_params()
                    return True
                else:
                    return False  # User cancelled the dialog
            except Exception as e:
                self.signalError.emit(f"Error reading NPZ file: {str(e)}")
                return False
        return True

    def browse_datafile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Data")
        if f:
            self.datafile_edit.setText(f)

    def browse_lppfile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select LPP")
        if f:
            self.filep_edit.setText(f)
