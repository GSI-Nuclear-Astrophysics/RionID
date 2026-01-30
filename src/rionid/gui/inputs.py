from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, 
                             QHBoxLayout, QComboBox, QCheckBox, QFileDialog, 
                             QMessageBox, QGroupBox, QApplication, QScrollArea, QToolButton, QFormLayout)
from PyQt5.QtCore import pyqtSignal, Qt, pyqtSlot
from PyQt5.QtGui import QFont, QCursor
import toml
import argparse
import logging as log
import numpy as np
import os
import time
import sys

from rionid.core import ImportData
from .controller import import_controller
from .dialogs import KeySelectionDialog

log.basicConfig(level=log.DEBUG)

class RionID_GUI(QWidget):
    """
    The main input control panel for RionID.

    This widget handles file selection, parameter configuration, and the execution
    of simulation scripts (Single Run and Quick PID). It communicates with the
    visualization widget to handle cursor picking and plot updates.
    """
    
    visualization_signal = pyqtSignal(object)
    overlay_sim_signal = pyqtSignal(object)
    clear_sim_signal = pyqtSignal()
    signalError = pyqtSignal(str)

    def __init__(self, plot_widget=None):
        super().__init__()
        self.visualization_widget = plot_widget
        self._stop_quick_pid = False
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
        self.setup_quick_pid()
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
                self.matched_result_edit.setText(p.get('matched_result', ''))
                self.alphap_min_edit.setText(p.get('alphap_min', ''))
                self.alphap_max_edit.setText(p.get('alphap_max', ''))
                self.alphap_step_edit.setText(p.get('alphap_step', ''))
                self.fref_min_edit.setText(p.get('fref_min', ''))
                self.fref_max_edit.setText(p.get('fref_max', ''))
                self.threshold_edit.setText(p.get('threshold', '1000'))
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
            'simulation_result': self.simulation_result_edit.text(),
            'matched_result': self.matched_result_edit.text(),
            'alphap_min': self.alphap_min_edit.text(),
            'alphap_max': self.alphap_max_edit.text(),
            'alphap_step': self.alphap_step_edit.text(),
            'fref_min': self.fref_min_edit.text(),
            'fref_max': self.fref_max_edit.text(),
            'threshold': self.threshold_edit.text()
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
        
        # Threshold for PID
        self.threshold_edit = QLineEdit("1000")
        h_t = QHBoxLayout()
        h_t.addWidget(QLabel("Match Threshold (Hz):"))
        h_t.addWidget(self.threshold_edit)
        self.vbox.addLayout(h_t)
        
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
        self.matched_result_edit = QLineEdit()
        opt_layout.addRow("Sim Result File:", self.simulation_result_edit)
        opt_layout.addRow("Matched Result File:", self.matched_result_edit)
        
        self.optional_group.setLayout(opt_layout)
        self.vbox.addWidget(self.optional_group)

    def setup_quick_pid(self):
        group = QGroupBox("Quick PID")
        layout = QVBoxLayout()
        
        self.alphap_min_edit = QLineEdit()
        self.alphap_max_edit = QLineEdit()
        self.alphap_step_edit = QLineEdit()
        
        h_a = QHBoxLayout()
        h_a.addWidget(QLabel("Alpha Range:"))
        h_a.addWidget(self.alphap_min_edit)
        h_a.addWidget(self.alphap_max_edit)
        h_a.addWidget(self.alphap_step_edit)
        layout.addLayout(h_a)
        
        self.fref_min_edit = QLineEdit()
        self.fref_max_edit = QLineEdit()
        
        pick_min = QPushButton("Pick")
        pick_min.clicked.connect(lambda: self.enterPlotPickMode(self.fref_min_edit))
        pick_max = QPushButton("Pick")
        pick_max.clicked.connect(lambda: self.enterPlotPickMode(self.fref_max_edit))
        
        h_f = QHBoxLayout()
        h_f.addWidget(QLabel("Freq Range:"))
        h_f.addWidget(self.fref_min_edit)
        h_f.addWidget(pick_min)
        h_f.addWidget(self.fref_max_edit)
        h_f.addWidget(pick_max)
        layout.addLayout(h_f)
        
        btn_pid = QPushButton("Run Quick PID")
        btn_pid.clicked.connect(self.quick_pid_script)
        layout.addWidget(btn_pid)
        
        group.setLayout(layout)
        self.vbox.addWidget(group)

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

    @pyqtSlot()
    def onPlotClicked(self):
        self._stop_quick_pid = True

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

    def quick_pid_script(self):
        """
        Executes the iterative Quick PID scanning algorithm.

        This method scans a range of Alpha_p values and Reference Frequencies (derived
        from experimental peaks) to find the best match (lowest Chi-squared) between
        simulation and experiment.
        """
        try:
            print("Running optimized quick_pid_script...")
            datafile = self.datafile_edit.text().strip()
            if not datafile: raise ValueError("No experimental data provided.")
            
            if not self._check_io_params(datafile): return

            # --- 1. Gather Constants & Ranges ---
            filep = self.filep_edit.text() or None
            remove_baseline = self.remove_baseline_checkbox.isChecked()
            psd_l = self._get_float(self.psd_baseline_removed_l_edit, 1000000.0)
            peak_pct = self._get_float(self.peak_thresh_edit, 0.05)
            min_dist = self._get_float(self.min_distance_edit, 10.0)
            harmonics = self.harmonics_edit.text()
            refion = self.refion_edit.text()
            circ = self._get_float(self.circumference_edit, 0.0)
            sim_sf = float(self.sim_scalingfactor_edit.text()) if self.sim_scalingfactor_edit.text() else None
            reload = self.reload_data_checkbox.isChecked()
            
            # Ranges
            f_min = self._get_float(self.fref_min_edit, -float('inf'))
            f_max = self._get_float(self.fref_max_edit, float('inf'))
            alpha_min = self._get_float(self.alphap_min_edit)
            alpha_max = self._get_float(self.alphap_max_edit)
            alpha_step = self._get_float(self.alphap_step_edit)
            threshold = self._get_float(self.threshold_edit, 1000.0)

            if alpha_step <= 0: alpha_step = 1e-5

            # --- 2. HEAVY LIFTING (Done ONCE) ---
            print("Loading data and calculating baseline (once)...")
            
            model = ImportData(
                refion=refion, alphap=alpha_min, filename=datafile, 
                circumference=circ, remove_baseline=remove_baseline, 
                psd_baseline_removed_l=psd_l, peak_threshold_pct=peak_pct, 
                min_distance=min_dist, io_params=self.current_io_params,
                reload_data=reload
            )
            
            model._set_particles_to_simulate_from_file(filep)
            model._calculate_moqs()
            
            if not hasattr(model, 'peak_freqs') or len(model.peak_freqs) == 0:
                raise RuntimeError("No experimental peaks detected.")
            
            exp_peaks = [f for f in model.peak_freqs if f_min <= f <= f_max]
            if not exp_peaks:
                raise RuntimeError("No peaks found in the specified Freq Range.")

            if isinstance(harmonics, str):
                harm_list = [float(h.strip()) for h in harmonics.replace(',', ' ').split()]
            else:
                harm_list = [float(harmonics)]

            # --- 3. THE FAST LOOP ---
            print(f"Scanning {len(exp_peaks)} peaks x Alpha range...")
            self._stop_quick_pid = False
            results = []
            
            orig_val_style = self.value_edit.styleSheet()
            orig_alpha_style = self.alphap_edit.styleSheet()

            for f_ref in exp_peaks:
                if self._stop_quick_pid: break
                self.value_edit.setText(f"{f_ref:.2f}")
                self.value_edit.setStyleSheet("background-color: #fff8b0;")
                QApplication.processEvents()

                for alpha in np.arange(alpha_min, alpha_max + 1e-12, alpha_step):
                    if self._stop_quick_pid: break
                    
                    # Update model parameters directly
                    model.alphap = alpha
                    # Recalculate Physics
                    model._calculate_srrf(fref=f_ref)
                    model._simulated_data(harmonics=harm_list, mode='Frequency', sim_scalingfactor=sim_sf)
                    
                    # Compute Match
                    chi2, count, _ = model.compute_matches(threshold)
                    results.append((f_ref, alpha, chi2, count))

            self.value_edit.setStyleSheet(orig_val_style)
            self.alphap_edit.setStyleSheet(orig_alpha_style)

            # --- 4. APPLY BEST RESULT ---
            if results:
                # Sort by Count (desc), then Chi2 (asc)
                best = sorted(results, key=lambda x: (-x[3], x[2]))[0]
                best_f, best_a, best_chi, best_cnt = best
                
                print(f"Best: F={best_f:.2f}, A={best_a:.6f}, Matches={best_cnt}, Chi2={best_chi:.2f}")
                self.value_edit.setText(f"{best_f:.2f}")
                self.alphap_edit.setText(f"{best_a:.6f}")
                
                model.alphap = best_a
                model._calculate_srrf(fref=best_f)
                model._simulated_data(harmonics=harm_list, mode='Frequency', sim_scalingfactor=sim_sf)
                model.compute_matches(threshold)
                
                self.saved_data = model
                self.visualization_signal.emit(model)
                self.save_parameters()

        except Exception as e:
            self.signalError.emit(str(e))
            import traceback
            traceback.print_exc()

    def browse_datafile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Data")
        if f: self.datafile_edit.setText(f)
    def browse_lppfile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select LPP")
        if f: self.filep_edit.setText(f)

class CollapsibleGroupBox(QGroupBox):
        def __init__(self, title="", parent=None):
            super().__init__(parent)
            self.setTitle("")
            self.toggle_button = QToolButton(text=title, checkable=True, checked=False)
            self.toggle_button.setStyleSheet("QToolButton { border: none; }")
            self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
            self.toggle_button.setArrowType(Qt.RightArrow)
            self.toggle_button.pressed.connect(self.on_pressed)

            self.content_widget = QWidget()
            self.content_layout = QVBoxLayout()
            self.content_layout.setContentsMargins(0, 0, 0, 0)
            self.content_widget.setLayout(self.content_layout)
            self.content_widget.setVisible(False)

            main_layout = QVBoxLayout()
            main_layout.addWidget(self.toggle_button)
            main_layout.addWidget(self.content_widget)
            main_layout.setContentsMargins(0, 0, 0, 0)
            self.setLayout(main_layout)

        def on_pressed(self):
            if self.toggle_button.isChecked():
                self.toggle_button.setArrowType(Qt.DownArrow)
                self.content_widget.setVisible(True)
            else:
                self.toggle_button.setArrowType(Qt.RightArrow)
                self.content_widget.setVisible(False)