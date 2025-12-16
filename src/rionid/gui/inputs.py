from PyQt5.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QCheckBox, QFileDialog, QMessageBox, QGroupBox, QToolButton
from PyQt5.QtCore import pyqtSignal, QThread, Qt
from PyQt5.QtGui import QFont
import toml
import argparse
import logging as log
from loguru import logger
from rionidgui.gui_controller import import_controller
import sys

log.basicConfig(level=log.DEBUG)
common_font = QFont()
common_font.setPointSize(12) #font size

class RionID_GUI(QWidget):
    visualization_signal = pyqtSignal(object)

    def __init__(self):
        super().__init__()
        self.initUI()
        self.load_parameters()  # Load parameters after initializing UI

    def initUI(self):
        self.setup_layout()
        self.setLayout(self.vbox)

    def load_parameters(self, filepath='parameters_cache.toml'):
        try:
            with open(filepath, 'r') as f:
                parameters = toml.load(f)
                self.datafile_edit.setText(parameters.get('datafile', ''))
                self.filep_edit.setText(parameters.get('filep', ''))
                self.alphap_edit.setText(parameters.get('alphap', ''))
                self.harmonics_edit.setText(parameters.get('harmonics', ''))
                self.refion_edit.setText(parameters.get('refion', ''))
                self.circumference_edit.setText(parameters.get('circumference', ''))
                self.mode_combo.setCurrentText(parameters.get('mode', 'Frequency'))
                self.value_edit.setText(parameters.get('value', ''))
                self.reload_data_checkbox.setChecked(parameters.get('reload_data', True))
                self.nions_edit.setText(parameters.get('nions', ''))
        except FileNotFoundError:
            pass  # No parameters file exists yet

    def save_parameters(self, filepath='parameters_cache.toml'):
        parameters = {
            'datafile': self.datafile_edit.text(),
            'filep': self.filep_edit.text(),
            'alphap': self.alphap_edit.text(),
            'harmonics': self.harmonics_edit.text(),
            'refion': self.refion_edit.text(),
            'circumference': self.circumference_edit.text(),
            'mode': self.mode_combo.currentText(),
            'value': self.value_edit.text(),
            'reload_data': self.reload_data_checkbox.isChecked(),
            'nions': self.nions_edit.text()
        }
        with open(filepath, 'w') as f:
            toml.dump(parameters, f)

    def setup_layout(self):
        self.vbox = QVBoxLayout()
        self.setup_file_selection()
        self.setup_parameters()
        self.setup_controls()

    def setup_file_selection(self):
        self.datafile_label = QLabel('Experimental Data File:')
        self.datafile_edit = QLineEdit()
        self.datafile_label.setFont(common_font)
        self.datafile_edit.setFont(common_font)
        self.datafile_button = QPushButton('Browse')
        self.datafile_button.setFont(common_font)

        self.datafile_button.clicked.connect(self.browse_datafile)

        self.filep_label = QLabel('.lpp File:')
        self.filep_edit = QLineEdit()
        self.filep_label.setFont(common_font)
        self.filep_edit.setFont(common_font)
        self.filep_button = QPushButton('Browse')
        self.filep_button.setFont(common_font)

        self.filep_button.clicked.connect(self.browse_lppfile)

        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.datafile_label)
        hbox1.addWidget(self.datafile_edit)
        hbox1.addWidget(self.datafile_button)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.filep_label)
        hbox2.addWidget(self.filep_edit)
        hbox2.addWidget(self.filep_button)

        self.vbox.addLayout(hbox1)
        self.vbox.addLayout(hbox2)

    def setup_parameters(self):
        self.alphap_label = QLabel('<i>&alpha;<sub>p</sub> or &gamma;<sub>t</sub> :</i>')
        self.alphap_edit = QLineEdit()
        self.alphap_label.setFont(common_font)
        self.alphap_edit.setFont(common_font)
        
        self.harmonics_label = QLabel('Harmonics (e.g 124 125 126):')
        self.harmonics_edit = QLineEdit()
        self.harmonics_label.setFont(common_font)
        self.harmonics_edit.setFont(common_font)

        self.refion_label = QLabel('Reference ion with format AAEl+QQ (e.g. 72Ge+32):')
        self.refion_edit = QLineEdit()
        self.refion_label.setFont(common_font)
        self.refion_edit.setFont(common_font)

        self.mode_label = QLabel('Mode:')
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(['Frequency', 'Bρ', 'Kinetic Energy'])
        self.mode_label.setFont(common_font)
        self.mode_combo.setFont(common_font)

        self.value_label = QLabel('Value:')
        self.value_edit = QLineEdit()
        self.value_label.setFont(common_font)
        self.value_edit.setFont(common_font)
        
        self.reload_data_checkbox = QCheckBox('Reload Experimental Data')
        self.reload_data_checkbox.setFont(common_font)
        self.reload_data_checkbox.setChecked(True)  # Default is to reload data

        self.circumference_label = QLabel('Circumference (m):')  # Add label for orbit length
        self.circumference_edit = QLineEdit()  # Add QLineEdit for orbit length  
        self.circumference_label.setFont(common_font)
        self.circumference_edit.setFont(common_font)      
        
        self.nions_label = QLabel('Number of ions to display (e.g. if 5, it will display the 5 most expected fragments):')
        self.nions_edit = QLineEdit()
        self.nions_label.setFont(common_font)
        self.nions_edit.setFont(common_font)    
        self.correction_label = QLabel('Add second order correction to the simulated frequencies with the form: <i>a<sub>0</sub> &middot; x<sup>2</sup> + a<sub>1</sub> &middot; x<sup>1</sup> + a<sub>2</sub> &middot; x<sup>0</sup>:</i>')        
        self.correction_edit = QLineEdit()
        self.correction_label.setFont(common_font)
        self.correction_edit.setFont(common_font)  
        
        hbox_alphap = QHBoxLayout()
        hbox_alphap.addWidget(self.alphap_label)
        hbox_alphap.addWidget(self.alphap_edit)

        hbox_harmonics = QHBoxLayout()
        hbox_harmonics.addWidget(self.harmonics_label)
        hbox_harmonics.addWidget(self.harmonics_edit)
        
        hbox_refion = QHBoxLayout()
        hbox_refion.addWidget(self.refion_label)
        hbox_refion.addWidget(self.refion_edit)

        hbox_circumference = QHBoxLayout()  # Add layout for orbit length
        hbox_circumference.addWidget(self.circumference_label)
        hbox_circumference.addWidget(self.circumference_edit)
        
        hbox_mode_value = QHBoxLayout()
        hbox_mode_value.addWidget(self.mode_label)
        hbox_mode_value.addWidget(self.mode_combo)
        hbox_mode_value.addWidget(self.value_edit)
        
        self.vbox.addWidget(self.reload_data_checkbox)
        self.vbox.addLayout(hbox_circumference)  # Add orbit length layout to the main layout
        self.vbox.addLayout(hbox_alphap)
        self.vbox.addLayout(hbox_harmonics)
        self.vbox.addLayout(hbox_refion)
        self.vbox.addLayout(hbox_mode_value)
        
        self.optional_features_group = CollapsibleGroupBox("Optional Features")
        self.optional_features_group.setFont(common_font)
        self.optional_features_group.addWidget(self.nions_label)
        self.optional_features_group.addWidget(self.nions_edit)
        self.optional_features_group.addWidget(self.correction_label)
        self.optional_features_group.addWidget(self.correction_edit)
        
        self.vbox.addWidget(self.optional_features_group)

    def setup_controls(self):
        self.run_button = QPushButton('Run')
        self.run_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        self.run_button.clicked.connect(self.run_script)


        self.exit_button = QPushButton('Exit')
        self.exit_button.setStyleSheet("""
            QPushButton {
                background-color: #f44336;
                color: white;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #d32f2f;
            }
        """)
        self.exit_button.clicked.connect(self.close_application)
        
        hbox_buttons = QHBoxLayout()
        hbox_buttons.addWidget(self.run_button)
        hbox_buttons.addWidget(self.exit_button)
        self.vbox.addLayout(hbox_buttons)
        
    def close_application(self):
        sys.exit()
        
    def browse_datafile(self):
        options = QFileDialog.Options()
        datafile, _ = QFileDialog.getOpenFileName(self, "Select Data File", "", "All Files (*);;NPZ Files (*.npz)", options= options)
        if datafile:
            self.datafile_edit.setText(datafile)

    def browse_lppfile(self):
        options = QFileDialog.Options()
        lppfile, _ = QFileDialog.getOpenFileName(self, "Select .lpp File", "", "All Files (*);;LPP Files (*.lpp)", options= options)
        if lppfile:
            self.filep_edit.setText(lppfile)

    def run_script(self):
        try:
            print("Running script...")
            datafile = self.datafile_edit.text()
            if not datafile:
                raise ValueError("No experimental data provided. Please enter any filename and click Run, the program will automatically calculate the simulated data.")

            filep = self.filep_edit.text()
            alphap = float(self.alphap_edit.text())
            harmonics = self.harmonics_edit.text()
            refion = self.refion_edit.text()
            circumference = float(self.circumference_edit.text())
            mode = self.mode_combo.currentText()
            value = self.value_edit.text()
            reload_data = self.reload_data_checkbox.isChecked()
            nions = self.nions_edit.text()

            args = argparse.Namespace(datafile=datafile,
                                        filep=filep or None,
                                        alphap=alphap or None,
                                        harmonics=harmonics or None,
                                        refion=refion or None,
                                        nions=nions or None,
                                        circumference=circumference or None,
                                        mode=mode or None,
                                        value=value or None,
                                        reload_data=reload_data or None)
            self.save_parameters()  # Save parameters before running the script

            # Simulate controller execution and emit data
            data = import_controller(**vars(args))
            if data:
                self.visualization_signal.emit(data)        
    
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'An error occurred: {str(e)}')
            log.error("Processing failed", exc_info=True)
            self.signalError.emit(str(e))

class CollapsibleGroupBox(QGroupBox):
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
