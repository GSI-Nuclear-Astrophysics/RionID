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
        self.setGeometry(10, 10, int(width * 0.9), int(height * 0.6))

        # Create a QSplitter to hold both the input and the visualization
        #splitter = QSplitter(Qt.Horizontal)
        #self.visualization_widget = CreatePyGUI()  # Initially empty
        
        #self.rion_input = RionID_GUI(plot_widget=self.visualization_widget)
        # Left panel - Input widget (RionID_GUI content)
        #splitter.addWidget(self.rion_input)
        # Right panel - Visualization widget (CreatePyGUI content)
        
        #splitter.addWidget(self.visualization_widget)
        
        # clear out old simulated curves when RionID_GUI says so
        #self.rion_input.clear_sim_signal.connect(
        #    self.visualization_widget.clear_simulated_data
        #)
        
        # Set initial size ratios (% input, % visualization)
        #splitter.setSizes([int(0.1*width), int(0.9*width)])
        # Dynamically resize both widgets
        #splitter.setStretchFactor(0, 1)  
        #splitter.setStretchFactor(1, 2) 
        
        # Create the main layout
        #layout = QVBoxLayout()
        #layout.addWidget(splitter)
        #self.setLayout(layout)
        
        # connect plot‐click signal to RionID_GUI’s stop slot
        #self.visualization_widget.plotClicked.connect(self.rion_input.onPlotClicked)
        
        # Connect the RionID_GUI signal to update CreatePyGUI once data is available
        #self.rion_input.visualization_signal.connect(self.update_visualization)
        # Connect the “overlay one simulation” signal so we don’t clear between curves
        #self.rion_input.overlay_sim_signal.connect(self.overlay_simulation)
        
    def overlay_simulation(self, data):
        """
        Receives one ImportData and overlays its simulated curves
        on the existing plot.  Clears out the old simulation items
        first so only the current batch is visible.
        """
        # 1) wipe out any previous simulated lines & labels
        self.visualization_widget.clear_simulated_data()
        # 2) draw the new ones
        self.visualization_widget.plot_simulated_data(data)
        # 3) force a repaint so we actually see it immediately
        from PyQt5.QtWidgets import QApplication
        QApplication.processEvents()   
      
    def update_visualization(self, data):
        """This method updates the visualization panel with new data."""
        self.visualization_widget.updateData(data)
