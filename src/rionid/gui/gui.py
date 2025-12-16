import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QMessageBox, QComboBox, QGroupBox, QGridLayout, QDesktopWidget,  QCheckBox, QSplitter
from PyQt5.QtCore import Qt, QLoggingCategory, QThread, pyqtSignal, QTimer, QEvent
import argparse
import os
import logging as log
from loguru import logger
from numpy import argsort, where, append, shape 
from rionid.importdata import ImportData
from rionid.pyqtgraphgui import CreatePyGUI
import toml
from barion.amedata import AMEData
from rionidgui.parameter_gui import RionID_GUI

log.basicConfig(level=log.DEBUG)

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RionID+")
        width = QDesktopWidget().screenGeometry(-1).width()
        height = QDesktopWidget().screenGeometry(-1).height()
        self.setGeometry(100, 100, width, height)  # Set window size

        # Create a QSplitter to hold both the input and the visualization
        splitter = QSplitter(Qt.Horizontal)

        # Left panel - Input widget (RionID_GUI content)
        self.rion_input = RionID_GUI()
        splitter.addWidget(self.rion_input)

        # Right panel - Visualization widget (CreatePyGUI content)
        self.visualization_widget = CreatePyGUI()  # Initially empty
        splitter.addWidget(self.visualization_widget)

        # Set initial size ratios (% input, % visualization)
        #splitter.setSizes([int(0.1*width), int(0.9*width)])
        # Dynamically resize both widgets
        splitter.setStretchFactor(0, 1)  
        splitter.setStretchFactor(1, 2) 

        # Create the main layout
        layout = QVBoxLayout()
        layout.addWidget(splitter)
        self.setLayout(layout)

        # Connect the RionID_GUI signal to update CreatePyGUI once data is available
        self.rion_input.visualization_signal.connect(self.update_visualization)

    def update_visualization(self, data):
        """This method updates the visualization panel with new data."""
        self.visualization_widget.updateData(data)