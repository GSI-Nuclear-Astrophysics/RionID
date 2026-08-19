from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout,
                             QHBoxLayout, QComboBox, QCheckBox, QFileDialog,
                             QMessageBox, QGroupBox, QScrollArea, QFormLayout)
from PyQt5.QtCore import pyqtSignal, Qt, pyqtSlot
from PyQt5.QtGui import QFont, QCursor
import toml
import argparse
import logging as log
import numpy as np
import os
import time
import sys

from .controller import import_controller
from .dialogs import KeySelectionDialog

log.basicConfig(level=log.DEBUG)

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

    def load_parameters(self, filepath='parameters_cache.toml'):
        try:
            with open(filepath, 'r') as f:
                p = toml.load(f)
                self.datafile_edit.setText(p.get('datafile', ''))
                self.filep_edit.setText(p.get('filep', ''))
                self.alphap_edit.setText(p.get('alphap', ''))
                self.harmonics_edit.setText(p.get('harmonics', ''))
                self.refion_edit.setText(p.get('refion', ''))
                self.circumference_edit.setText(p.get('circumference', ''))
                self.highlight_ions_edit.setText(p.get('highlight_ions', ''))
                self.mode_combo.setCurrentText(p.get('mode', 'Frequency'))
                self.value_edit.setText(p.get('value', ''))
                self.sim_scalingfactor_edit.setText(p.get('sim_scalingfactor', ''))
                self.remove_baseline_checkbox.setChecked(p.get('remove_baseline_checkbox', False))
                self.psd_baseline_removed_l_edit.setText(str(p.get('psd_baseline_removed_l', '1000000')))
                self.peak_thresh_edit.setText(str(p.get('peak_threshold_pct', '0.05')))
                self.min_distance_edit.setText(str(p.get('min_distance', '10')))
                self.matching_freq_min_edit.setText(str(p.get('matching_freq_min', '')))
                self.matching_freq_max_edit.setText(str(p.get('matching_freq_max', '')))
                self.correction_edit.setText(p.get('correction', ''))
                self.nions_edit.setText(p.get('nions', ''))
                self.reload_data_checkbox.setChecked(p.get('reload_data', True))
                self.simulation_result_edit.setText(p.get('simulation_result', ''))
        except FileNotFoundError: pass

    def save_parameters(self, filepath='parameters_cache.toml'):
        p = {
            'datafile': self.datafile_edit.text(),
            'filep': self.filep_edit.text(),
            'alphap': self.alphap_edit.text(),
            'harmonics': self.harmonics_edit.text(),
            'refion': self.refion_edit.text(),
            'circumference': self.circumference_edit.text(),
            'highlight_ions': self.highlight_ions_edit.text(),
            'mode': self.mode_combo.currentText(),
            'value': self.value_edit.text(),
            'sim_scalingfactor': self.sim_scalingfactor_edit.text(),
            'remove_baseline_checkbox': self.remove_baseline_checkbox.isChecked(),
            'psd_baseline_removed_l': self.psd_baseline_removed_l_edit.text(),
            'peak_threshold_pct': self.peak_thresh_edit.text(),
            'min_distance': self.min_distance_edit.text(),
            'matching_freq_min': self.matching_freq_min_edit.text(),
            'matching_freq_max': self.matching_freq_max_edit.text(),
            'correction': self.correction_edit.text(),
            'nions': self.nions_edit.text(),
            'reload_data': self.reload_data_checkbox.isChecked(),
            'simulation_result': self.simulation_result_edit.text()
        }
        with open(filepath, 'w') as f: toml.dump(p, f)

    def setup_file_selection(self):
        self.datafile_label = QLabel('Experimental Data File:')
        self.datafile_edit = QLineEdit()
        self.datafile_button = QPushButton('Browse')
        self.datafile_button.clicked.connect(self.browse_datafile)
        
        self.filep_label = QLabel('.lpp File:')
        self.filep_edit = QLineEdit()
        self.filep_button = QPushButton('Browse')
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
        # Baseline
        self.remove_baseline_checkbox = QCheckBox('Remove Baseline')
        self.vbox.addWidget(self.remove_baseline_checkbox)
        
        self.psd_baseline_removed_l_edit = QLineEdit("1000000")
        hb_bl = QHBoxLayout()
        hb_bl.addWidget(QLabel("Baseline l:"))
        hb_bl.addWidget(self.psd_baseline_removed_l_edit)
        self.vbox.addLayout(hb_bl)

        # Alpha P
        self.alphap_edit = QLineEdit()
        hb_ap = QHBoxLayout()
        hb_ap.addWidget(QLabel("Alpha P:"))
        hb_ap.addWidget(self.alphap_edit)
        self.vbox.addLayout(hb_ap)

        # Standard Params
        self.harmonics_edit = QLineEdit()
        self.refion_edit = QLineEdit()
        self.circumference_edit = QLineEdit()
        self.highlight_ions_edit = QLineEdit()
        
        for lbl, widget in [("Harmonics:", self.harmonics_edit), 
                            ("Ref Ion:", self.refion_edit),
                            ("Circumference:", self.circumference_edit),
                            ("Highlight Ions:", self.highlight_ions_edit)]:
            h = QHBoxLayout()
            h.addWidget(QLabel(lbl))
            h.addWidget(widget)
            self.vbox.addLayout(h)

        # Mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(['Frequency', 'Bρ', 'Kinetic Energy'])
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

        # Peak Detection
        self.peak_thresh_edit = QLineEdit("0.05")
        self.min_distance_edit = QLineEdit("10")
        h_peak = QHBoxLayout()
        h_peak.addWidget(QLabel("Peak Thresh %:"))
        h_peak.addWidget(self.peak_thresh_edit)
        h_peak.addWidget(QLabel("Min Dist:"))
        h_peak.addWidget(self.min_distance_edit)
        self.vbox.addLayout(h_peak)
        
        # Matching Freq Range
        self.matching_freq_min_edit = QLineEdit()
        self.matching_freq_max_edit = QLineEdit()
        self.pick_matching_freq_min_button = QPushButton("Pick")
        self.pick_matching_freq_min_button.clicked.connect(lambda: self.enterPlotPickMode(self.matching_freq_min_edit))
        self.pick_matching_freq_max_button = QPushButton("Pick")
        self.pick_matching_freq_max_button.clicked.connect(lambda: self.enterPlotPickMode(self.matching_freq_max_edit))
        
        h_mf = QHBoxLayout()
        h_mf.addWidget(QLabel("Match Freq Min:"))
        h_mf.addWidget(self.matching_freq_min_edit)
        h_mf.addWidget(self.pick_matching_freq_min_button)
        self.vbox.addLayout(h_mf)
        
        h_mf2 = QHBoxLayout()
        h_mf2.addWidget(QLabel("Match Freq Max:"))
        h_mf2.addWidget(self.matching_freq_max_edit)
        h_mf2.addWidget(self.pick_matching_freq_max_button)
        self.vbox.addLayout(h_mf2)
        
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
        self.run_button.clicked.connect(self.run_script)
        self.vbox.addWidget(self.run_button)
        
        self.exit_button = QPushButton("Exit")
        self.exit_button.clicked.connect(self.close_application)
        self.vbox.addWidget(self.exit_button)

    def close_application(self):
        sys.exit()

    def enterPlotPickMode(self, target):
        if not self.visualization_widget: return
        self._pick_target = target
        target.setStyleSheet("background-color: lightgray;")
        self.visualization_widget.plot_widget.setCursor(Qt.CrossCursor)
        self.visualization_widget.plotClicked.connect(self._onPlotPicked)

    @pyqtSlot()
    def _onPlotPicked(self):
        pos = self.visualization_widget.plot_widget.mapFromGlobal(QCursor.pos())
        point = self.visualization_widget.plot_widget.plotItem.vb.mapSceneToView(pos)
        
        if self._pick_target:
            self._pick_target.setText(f"{point.x()*1e6:.2f}") 
            self._pick_target.setStyleSheet("")
        
        self.visualization_widget.plot_widget.setCursor(Qt.ArrowCursor)
        self.visualization_widget.plotClicked.disconnect(self._onPlotPicked)
    
    def _get_float(self, widget, default=0.0):
        text = widget.text().strip()
        if not text: return default
        try: return float(text)
        except ValueError: return default

    def run_script(self):
        datafile = self.datafile_edit.text()
        if not datafile: return
        
        io_params = self.current_io_params
        ext = os.path.splitext(datafile)[1].lower()
        if ext == '.npz' and not io_params:
            data = np.load(datafile)
            dlg = KeySelectionDialog(self, list(data.keys()))
            if dlg.exec_(): 
                io_params = dlg.get_params()
                self.current_io_params = io_params 
            else: return

        correct_str = self.correction_edit.text().strip()
        correct = [float(x) for x in correct_str.split()] if correct_str else None
        
        sim_sf_str = self.sim_scalingfactor_edit.text().strip()
        sim_sf = float(sim_sf_str) if sim_sf_str else None

        psd_l = self._get_float(self.psd_baseline_removed_l_edit, 1000000.0)
        peak_pct = self._get_float(self.peak_thresh_edit, 0.05)
        min_dist = self._get_float(self.min_distance_edit, 10.0)
        alphap = self._get_float(self.alphap_edit, 0.0)
        circumference = self._get_float(self.circumference_edit, 0.0)
        
        match_min_str = self.matching_freq_min_edit.text().strip()
        match_min = float(match_min_str) if match_min_str else None
        match_max_str = self.matching_freq_max_edit.text().strip()
        match_max = float(match_max_str) if match_max_str else None

        args = argparse.Namespace(
            datafile=datafile,
            filep=self.filep_edit.text(),
            alphap=alphap,
            harmonics=self.harmonics_edit.text(),
            refion=self.refion_edit.text(),
            circumference=circumference,
            mode=self.mode_combo.currentText(),
            value=self.value_edit.text(),
            remove_baseline=self.remove_baseline_checkbox.isChecked(),
            psd_baseline_removed_l=psd_l,
            peak_threshold_pct=peak_pct,
            min_distance=min_dist,
            highlight_ions=self.highlight_ions_edit.text(),
            io_params=io_params,
            reload_data=self.reload_data_checkbox.isChecked(),
            nions=self.nions_edit.text(),
            sim_scalingfactor=sim_sf,
            matching_freq_min=match_min,
            matching_freq_max=match_max,
            correct=correct
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
        if ext == '.npz' and not self.current_io_params:
            try:
                data = np.load(datafile)
                dlg = KeySelectionDialog(self, list(data.keys()))
                if dlg.exec_(): 
                    self.current_io_params = dlg.get_params()
                    return True
                else:
                    return False # User cancelled the dialog
            except Exception as e:
                self.signalError.emit(f"Error reading NPZ file: {str(e)}")
                return False
        return True
        
    def browse_datafile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Data")
        if f: self.datafile_edit.setText(f)
    def browse_lppfile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select LPP")
        if f: self.filep_edit.setText(f)