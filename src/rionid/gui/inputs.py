from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, 
                             QHBoxLayout, QComboBox, QCheckBox, QFileDialog, 
                             QMessageBox, QGroupBox, QApplication, QScrollArea, QToolButton)
from PyQt5.QtCore import pyqtSignal, Qt, pyqtSlot
from PyQt5.QtGui import QFont, QCursor
import toml
import argparse
import logging as log
import numpy as np
import os
import time

from rionid.core import ImportData
from .controller import import_controller
from .dialogs import KeySelectionDialog

log.basicConfig(level=log.DEBUG)
common_font = QFont()
common_font.setPointSize(12) 

class RionID_GUI(QWidget):
    """
    The main input control panel for RionID.

    This widget handles file selection, parameter configuration, and the execution
    of simulation scripts (Single Run and Quick PID). It communicates with the
    visualization widget to handle cursor picking and plot updates.

    Attributes
    ----------
    visualization_signal : pyqtSignal
        Emits the resulting `ImportData` object after a successful simulation.
    overlay_sim_signal : pyqtSignal
        Emits an `ImportData` object to overlay on the existing plot (used in PID loops).
    clear_sim_signal : pyqtSignal
        Signal to clear simulated lines from the plot.
    signalError : pyqtSignal
        Emits error messages to be displayed in a message box.
    """
    
    visualization_signal = pyqtSignal(object)
    overlay_sim_signal = pyqtSignal(object)
    clear_sim_signal = pyqtSignal()
    signalError = pyqtSignal(str)

    def __init__(self, plot_widget=None):
        """
        Initialize the GUI inputs panel.

        Parameters
        ----------
        plot_widget : CreatePyGUI, optional
            Reference to the plotting widget. Required for 'Pick' cursor functionality.
        """
        super().__init__()
        self.visualization_widget = plot_widget
        self._stop_quick_pid = False
        self.saved_data = None
        self.current_io_params = {} 
        self.initUI()
        self.load_parameters()
        self.signalError.connect(self.show_error)

    def show_error(self, msg):
        """Display a critical error message box."""
        QMessageBox.critical(self, "Error", msg)

    def initUI(self):
        """Sets up the scrollable layout and initializes all sub-components."""
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
        """
        Loads configuration parameters from a TOML file.

        Parameters
        ----------
        filepath : str
            Path to the cache file. Defaults to 'parameters_cache.toml'.
        """
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
        except FileNotFoundError:
            pass 
        except Exception as e:
            print(f"Error loading parameters: {e}")

    def save_parameters(self, filepath='parameters_cache.toml'):
        """
        Saves current GUI state to a TOML file.

        Parameters
        ----------
        filepath : str
            Path to the cache file.
        """
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
        with open(filepath, 'w') as f:
            toml.dump(p, f)

    def setup_file_selection(self):
        """Creates the file selection widgets (Data file, LISE++ file)."""
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
        """Creates the main physics parameter widgets."""
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
        self.optional_group = CollapsibleGroupBox("Optional Features")
        opt_layout = QVBoxLayout()
        
        self.nions_edit = QLineEdit()
        h_n = QHBoxLayout()
        h_n.addWidget(QLabel("N Ions to Display:"))
        h_n.addWidget(self.nions_edit)
        opt_layout.addLayout(h_n)
        
        self.correction_edit = QLineEdit()
        self.correction_edit.setPlaceholderText("a0 a1 a2")
        h_c = QHBoxLayout()
        h_c.addWidget(QLabel("Correction (a0 a1 a2):"))
        h_c.addWidget(self.correction_edit)
        opt_layout.addLayout(h_c)
        
        self.reload_data_checkbox = QCheckBox("Reload Data Cache")
        opt_layout.addWidget(self.reload_data_checkbox)
        
        self.simulation_result_edit = QLineEdit()
        self.matched_result_edit = QLineEdit()
        opt_layout.addWidget(QLabel("Sim Result File:"))
        opt_layout.addWidget(self.simulation_result_edit)
        opt_layout.addWidget(QLabel("Matched Result File:"))
        opt_layout.addWidget(self.matched_result_edit)
        
        self.optional_group.setLayout(opt_layout)
        self.vbox.addWidget(self.optional_group)

    def setup_quick_pid(self):
        """Creates the Quick PID scanning widgets."""
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
        """Creates the main Run/Exit buttons."""
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_script)
        self.vbox.addWidget(self.run_button)
        
        self.exit_button = QPushButton("Exit")
        self.exit_button.clicked.connect(self.close_application)
        self.vbox.addWidget(self.exit_button)

    def close_application(self):
        sys.exit()

    def enterPlotPickMode(self, target):
        """
        Enters 'Pick Mode' where the next click on the plot captures the X-coordinate.

        Parameters
        ----------
        target : QLineEdit
            The text field where the picked value will be inserted.
        """
        if not self.visualization_widget: return
        self._pick_target = target
        target.setStyleSheet("background-color: lightgray;")
        self.visualization_widget.plot_widget.setCursor(Qt.CrossCursor)
        self.visualization_widget.plotClicked.connect(self._onPlotPicked)

    @pyqtSlot()
    def _onPlotPicked(self):
        """Slot called when the plot is clicked in Pick Mode."""
        pos = self.visualization_widget.plot_widget.mapFromGlobal(QCursor.pos())
        point = self.visualization_widget.plot_widget.plotItem.vb.mapSceneToView(pos)
        
        if self._pick_target:
            self._pick_target.setText(f"{point.x()*1e6:.2f}") 
            self._pick_target.setStyleSheet("")
        
        self.visualization_widget.plot_widget.setCursor(Qt.ArrowCursor)
        self.visualization_widget.plotClicked.disconnect(self._onPlotPicked)

    @pyqtSlot()
    def onPlotClicked(self):
        """Slot called to interrupt the Quick PID loop."""
        self._stop_quick_pid = True

    def run_script(self):
        """
        Executes a single simulation run based on current parameters.
        """
        datafile = self.datafile_edit.text()
        if not datafile: return
        
        # Handle NPZ key selection
        io_params = self.current_io_params
        ext = os.path.splitext(datafile)[1].lower()
        if ext == '.npz' and not io_params:
            data = np.load(datafile)
            dlg = KeySelectionDialog(self, list(data.keys()))
            if dlg.exec_(): 
                io_params = dlg.get_params()
                self.current_io_params = io_params 
            else: return

        # Parse Correction
        correct_str = self.correction_edit.text().strip()
        correct = [float(x) for x in correct_str.split()] if correct_str else None
        
        # Parse Scaling
        sim_sf_str = self.sim_scalingfactor_edit.text().strip()
        sim_sf = float(sim_sf_str) if sim_sf_str else None

        args = argparse.Namespace(
            datafile=datafile,
            filep=self.filep_edit.text(),
            alphap=float(self.alphap_edit.text() or 0),
            harmonics=self.harmonics_edit.text(),
            refion=self.refion_edit.text(),
            circumference=float(self.circumference_edit.text() or 0),
            mode=self.mode_combo.currentText(),
            value=self.value_edit.text(),
            remove_baseline=self.remove_baseline_checkbox.isChecked(),
            psd_baseline_removed_l=float(self.psd_baseline_removed_l_edit.text()),
            peak_threshold_pct=float(self.peak_thresh_edit.text()),
            min_distance=float(self.min_distance_edit.text()),
            highlight_ions=self.highlight_ions_edit.text(),
            io_params=io_params,
            reload_data=self.reload_data_checkbox.isChecked(),
            nions=self.nions_edit.text(),
            sim_scalingfactor=sim_sf,
            matching_freq_min=float(self.matching_freq_min_edit.text() or 0) if self.matching_freq_min_edit.text() else None,
            matching_freq_max=float(self.matching_freq_max_edit.text() or 0) if self.matching_freq_max_edit.text() else None,
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
            print("Running quick_pid_script…")
            datafile = self.datafile_edit.text().strip()
            if not datafile:
                raise ValueError("No experimental data provided.")

            # --- collect constant arguments once ---
            filep = self.filep_edit.text() or None
            remove_baseline = self.remove_baseline_checkbox.isChecked()
            psd_baseline_removed_l = float(self.psd_baseline_removed_l_edit.text())
            alphap = float(self.alphap_edit.text())
            peak_threshold_pct = float(self.peak_thresh_edit.text())
            min_distance = float(self.min_distance_edit.text())
            harmonics = self.harmonics_edit.text()
            refion = self.refion_edit.text()
            highlight_ions = self.highlight_ions_edit.text() or None
            nions = self.nions_edit.text() or None
            circumference = float(self.circumference_edit.text())
            sim_scalingfactor = self.sim_scalingfactor_edit.text().strip()
            sim_scalingfactor = float(sim_scalingfactor) if sim_scalingfactor else None
            reload_data = self.reload_data_checkbox.isChecked()
            simulation_result= self.simulation_result_edit.text()
            matched_result= self.matched_result_edit.text()

            try:
                threshold = float(self.threshold_edit.text())
            except ValueError:
                raise ValueError("Please enter a valid number for matching threshold")
            try:
                matching_freq_min = float(self.matching_freq_min_edit.text())
                matching_freq_max = float(self.matching_freq_max_edit.text())
            except ValueError:
                # Allow empty range for initial detection
                matching_freq_min = None
                matching_freq_max = None
                
            fref_min = float(self.fref_min_edit.text() or '-inf')
            fref_max = float(self.fref_max_edit.text() or 'inf')

            # --- 1) Load experimental data and detect peaks ---
            model = ImportData(
                refion=refion,
                highlight_ions=highlight_ions,
                remove_baseline = remove_baseline or None,
                psd_baseline_removed_l=psd_baseline_removed_l or None,
                alphap=alphap,
                filename=datafile,
                reload_data=reload_data,
                circumference=circumference,
                peak_threshold_pct=peak_threshold_pct,
                min_distance=min_distance,
                matching_freq_min=matching_freq_min,
                matching_freq_max=matching_freq_max,
                io_params=self.current_io_params
            )
            if not hasattr(model, 'peak_freqs') or len(model.peak_freqs) == 0:
                raise RuntimeError("Could not detect any experimental peaks.")
            self.visualization_signal.emit(model)
            
            # experimental peak frequencies (Hz)
            exp_peaks_hz = model.peak_freqs
            print(f"Detected {len(exp_peaks_hz)} experimental peaks.")

            # define your alphap scan range
            alphap_min  = float(self.alphap_min_edit.text())
            alphap_max  = float(self.alphap_max_edit.text())
            alphap_step = float(self.alphap_step_edit.text())
            
            self._stop_quick_pid = False
            QApplication.processEvents()

            results = []
 
            # Filter experimental peaks
            exp_peaks_hz_filtering = [f for f in model.peak_freqs if fref_min <= f <= fref_max]
            
            if not exp_peaks_hz_filtering:
                QMessageBox.critical(self, "Error", "No experimental peaks found within the specified frequency range.")
                self.fref_min_edit.setStyleSheet("background-color: red;")
                self.fref_max_edit.setStyleSheet("background-color: red;")
                return
            else:
                self.fref_min_edit.setStyleSheet("")
                self.fref_max_edit.setStyleSheet("")
                   
            # Grab and remember the original styles
            orig_value_style  = self.value_edit.styleSheet()
            orig_alpha_style  = self.alphap_edit.styleSheet()
            
            # Initialize first iteration flag
            first_iteration = True
            
            for f_ref in exp_peaks_hz_filtering:
                QApplication.processEvents()
                if self._stop_quick_pid:
                    print("Quick‐PID scan was stopped by user click.")
                    break

                # Highlight current f_ref in the UI
                self.value_edit.setStyleSheet("background-color: #fff8b0;")  
                self.value_edit.setText(f"{f_ref:.2f}")
                QApplication.processEvents()

                # Inner loop over a range of test_alphap values
                for test_alphap in np.arange(alphap_min, alphap_max + 1e-12, alphap_step):
                    if self._stop_quick_pid: break
                    
                    # Update UI to show which alphap is being tested
                    self.alphap_edit.setStyleSheet("background-color: #b0fff8;")
                    self.alphap_edit.setText(f"{test_alphap:.6f}")
                    QApplication.processEvents()

                    # Run simulation for this combination
                    sim_args = argparse.Namespace(
                        datafile=datafile,
                        filep=filep,
                        remove_baseline = remove_baseline,
                        psd_baseline_removed_l = psd_baseline_removed_l,
                        alphap=test_alphap,
                        harmonics=harmonics,
                        refion=refion,
                        highlight_ions=highlight_ions,
                        nions=nions,
                        circumference=circumference,
                        mode='Frequency',
                        sim_scalingfactor=sim_scalingfactor,
                        value=f_ref,
                        reload_data=reload_data if first_iteration else False, # Cache optimization
                        peak_threshold_pct=peak_threshold_pct,
                        min_distance=min_distance,
                        output_results=False,
                        matching_freq_min=matching_freq_min,
                        matching_freq_max=matching_freq_max,
                        simulation_result=simulation_result,
                        io_params=self.current_io_params,
                        correct=None
                    )
                    
                    data_i = import_controller(**vars(sim_args))
                    if data_i is None: continue
                    
                    chi2, match_count, highlights = data_i.compute_matches(threshold, matching_freq_min, matching_freq_max)
                    results.append((f_ref, test_alphap, chi2, match_count, highlights))
                    
                    if first_iteration:
                        self.saved_data = data_i
                        self.overlay_sim_signal.emit(self.saved_data)
                        first_iteration = False
                    
                    del data_i

                if not results: continue
                
                # Sort results: Max Matches (desc), Min Chi2 (asc)
                sorted_results = sorted(results, key=lambda x: (-x[3], x[2]))
                best_fref, best_alphap, best_chi2, best_match_count, best_match_ions = sorted_results[0]
                
                # Run simulation for best combination
                sim_args = argparse.Namespace(
                    datafile=datafile,
                    filep=filep,
                    remove_baseline = remove_baseline,
                    psd_baseline_removed_l = psd_baseline_removed_l,
                    alphap=best_alphap,
                    harmonics=harmonics,
                    refion=refion,
                    highlight_ions=highlight_ions,
                    nions=nions,
                    circumference=circumference,
                    mode='Frequency',
                    sim_scalingfactor=sim_scalingfactor,
                    value=best_fref,
                    reload_data=False,
                    peak_threshold_pct=peak_threshold_pct,
                    min_distance=min_distance,
                    output_results=True,
                    matching_freq_min=matching_freq_min,
                    matching_freq_max=matching_freq_max,
                    simulation_result=simulation_result,
                    io_params=self.current_io_params,
                    correct=None
                )
                best_data = import_controller(**vars(sim_args))
                best_chi2, best_match_count, best_match_ions = best_data.compute_matches(threshold, matching_freq_min, matching_freq_max)
                best_data.save_matched_result(matched_result)
                
                self.save_parameters()
                print(f"\n→ Best: f_ref={best_fref:.2f}Hz, alphap={best_alphap:.4f}, χ²={best_chi2:.3e}, matches={best_match_count}")
                
                self.mode_combo.setCurrentText('Frequency')
                self.value_edit.setText(f"{best_fref:.2f}")
                self.alphap_edit.setText(f"{best_alphap:.6f}")
                self.overlay_sim_signal.emit(best_data)
                QApplication.processEvents()
                
                self.alphap_edit.setStyleSheet(orig_alpha_style)
                
            self.value_edit.setStyleSheet(orig_value_style)
            self.save_parameters()

        except Exception as e:
            QMessageBox.critical(self, "Quick PID Error", str(e))
            log.error("quick_pid_script failed", exc_info=True)

    def browse_datafile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select Data")
        if f: self.datafile_edit.setText(f)
    def browse_lppfile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select LPP")
        if f: self.filep_edit.setText(f)

class CollapsibleGroupBox(QGroupBox):
    """
    A custom QGroupBox that can be collapsed/expanded by clicking its title.
    """
    def __init__(self, title="", parent=None):
        super(CollapsibleGroupBox, self).__init__(parent)
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

    def addWidget(self, widget):
        self.content_layout.addWidget(widget)