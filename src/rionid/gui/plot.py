import sys
import numpy as np
import re
import pyqtgraph as pg
from PyQt5.QtWidgets import (QMainWindow, QApplication, QVBoxLayout, QWidget, 
                             QPushButton, QHBoxLayout, QLabel, QDesktopWidget, QSpinBox)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QLoggingCategory, Qt, pyqtSignal

class CustomLegendItem(pg.LegendItem):
    """
    A subclass of pyqtgraph.LegendItem that supports dynamic font sizing.

    Standard pyqtgraph legends do not easily support font updates after initialization.
    This class overrides the item addition to apply a custom QFont.

    Parameters
    ----------
    font_size : int
        The initial point size of the font.
    """
    def __init__(self, font_size, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.font = QFont("Arial", font_size)

    def addItem(self, item, name):
        """
        Adds an item to the legend with the custom font applied to the label.
        """
        label = pg.LabelItem(text=name, justify='left')
        label.setFont(self.font)
        super().addItem(item, name)
    
    def updateFont(self, font_size):
        """
        Updates the font size for future items added to the legend.
        
        Note: Existing items are not automatically resized in this simple implementation;
        the legend is usually cleared and rebuilt during plot updates.
        """
        self.font.setPointSize(font_size)

class CreatePyGUI(QMainWindow):
    """
    The main visualization widget for RionID.

    This class handles the plotting of experimental spectra (blue lines),
    detected peaks (red triangles), and simulated ion revolution frequencies
    (dashed lines). It supports zooming, cursor tracking, and dynamic font sizing.

    Attributes
    ----------
    plotClicked : pyqtSignal
        Signal emitted when the user clicks on the plot area (used for 'Pick' mode).
    plot_widget : pg.PlotWidget
        The central pyqtgraph widget.
    """
    plotClicked = pyqtSignal()

    def __init__(self, exp_data=None, sim_data=None):
        """
        Initializes the plotting window.

        Parameters
        ----------
        exp_data : tuple, optional
            Initial experimental data (freq, amp).
        sim_data : dict, optional
            Initial simulation dictionary.
        """
        super().__init__()
        self.saved_x_range = None  
        self.simulated_items = []
        self.red_triangles = None
        self.exp_data_line = None
        self.font_size = 20
        
        self.setup_ui()
        
        # Connect scene click for "Pick Mode"
        self.plot_widget.scene().sigMouseClicked.connect(self.on_click)

    def on_click(self, event):
        """Emits the plotClicked signal when the mouse is clicked."""
        self.plotClicked.emit()

    def setup_ui(self):
        """Sets up the layout, plot widget styling, axes, and control buttons."""
        self.setWindowTitle('Schottky Signals Identifier')
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        main_layout = QVBoxLayout(self.main_widget)
        
        # Suppress annoying Qt warnings
        QLoggingCategory.setFilterRules('*.warning=false\n*.critical=false')
        
        # Configure PlotWidget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w') # White background
        self.plot_widget.plotItem.ctrl.logYCheck.setChecked(True) # Default to Log Y
        
        # Configure Legend
        self.legend = CustomLegendItem(self.font_size, offset=(-10, 10))
        self.legend.setParentItem(self.plot_widget.graphicsItem())
        self.legend.setBrush(pg.mkBrush('white'))
        self.legend.setLabelTextColor('black')
        
        # Style Axes
        self.plot_widget.getAxis('bottom').setPen(pg.mkPen('black'))
        self.plot_widget.getAxis('left').setPen(pg.mkPen('black'))
        self.plot_widget.getAxis('bottom').setTextPen('black')
        self.plot_widget.getAxis('left').setTextPen('black')
        
        main_layout.addWidget(self.plot_widget)
        
        # Cursor Label
        self.cursor_pos_label = QLabel(self)
        self.cursor_pos_label.setStyleSheet("color: black;")
        main_layout.addWidget(self.cursor_pos_label)
        self.proxy = pg.SignalProxy(self.plot_widget.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        
        self.add_buttons(main_layout)
        self.update_fonts(self.font_size)

    def plot_all_data(self, data):
        """
        Clears the plot and redraws both experimental and simulated data.

        Parameters
        ----------
        data : ImportData
            The data model object containing experimental arrays and simulation results.
        """
        # We clear specific items rather than the whole widget to preserve Legend/Axis config
        self.clear_experimental_data()
        self.clear_simulated_data()
        
        self.plot_experimental_data(data)
        self.plot_simulated_data(data)

    def plot_experimental_data(self, data):
        """
        Plots the experimental spectrum and detected peaks.

        Parameters
        ----------
        data : ImportData
            Must contain `experimental_data` (tuple) and optionally `peak_freqs`.
        """
        if data.experimental_data is None: return
        self.exp_data = data.experimental_data
        
        self.x_exp, self.z_exp = self.exp_data[0]*1e-6, self.exp_data[1]
        
        # Auto-range only on first load
        if self.saved_x_range is None:
            self.saved_x_range = (min(self.x_exp), max(self.x_exp))
            self.plot_widget.setXRange(*self.saved_x_range, padding=0.05)

        # Plot Spectrum (Blue Line)
        self.exp_data_line = self.plot_widget.plot(self.x_exp, self.z_exp, pen=pg.mkPen('blue', width=2))
        self.legend.addItem(self.exp_data_line, 'Experimental Data')
        
        # Plot Peaks (Red Triangles)
        if hasattr(data, 'peak_freqs') and len(data.peak_freqs) > 0:
            self.red_triangles = self.plot_widget.plot(
                data.peak_freqs * 1e-6, data.peak_heights,
                pen=None, symbol='t', symbolBrush='r', symbolSize=12
            )
            self.legend.addItem(self.red_triangles, 'Peaks')

    def plot_simulated_data(self, data):
        """
        Plots vertical lines for simulated ion frequencies.

        Applies color coding:
        - Green: Highlighted ions (matches).
        - Orange: Reference ion.
        - Varied Colors: Other harmonics.

        Parameters
        ----------
        data : ImportData
            Must contain `simulated_data_dict`, `ref_ion`, and `highlight_ions`.
        """
        self.simulated_data = data.simulated_data_dict
        refion = data.ref_ion
        highlights = data.highlight_ions or []
        
        for i, (harmonic, sdata) in enumerate(self.simulated_data.items()):
            color = pg.intColor(i, hues=len(self.simulated_data))
            for entry in sdata:
                freq = float(entry[0])*1e-6
                label = entry[2]
                yield_value = float(entry[1])
                
                # Determine Color
                if label in highlights: label_color = 'green'
                elif label == refion: label_color = 'orange'
                else: label_color = color
                
                # Format Label (Superscript)
                match = re.match(r'(\d+)([A-Za-z]+)(\d+)\+', label)
                if match:
                    mass, elem, charge = match.groups()
                    new_label = self.to_superscript(mass) + elem + self.to_superscript(charge) + '⁺'
                else: new_label = label
                
                # Plot Line
                line = self.plot_widget.plot([freq, freq], [1e-10, yield_value], pen=pg.mkPen(color=label_color, width=1, style=Qt.DashLine))
                
                # Plot Label
                text = pg.TextItem(text=new_label, color=label_color, anchor=(0, 0.5))
                text.setFont(QFont("Arial", self.font_size))
                text.setAngle(90)
                text.setPos(freq, yield_value * 1.2)
                self.plot_widget.addItem(text)
                
                self.simulated_items.append((line, text))
            
            self.legend.addItem(line, f'Harmonic {harmonic}')

    def to_superscript(self, s):
        """Converts a string of numbers to unicode superscripts."""
        supers = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}
        return ''.join(supers.get(c, c) for c in s)

    def update_fonts(self, size):
        """Updates the font size for axes and legends dynamically."""
        self.font_size = size
        self.font_ticks = QFont()
        self.font_ticks.setPixelSize(size)
        self.plot_widget.getAxis('bottom').setTickFont(self.font_ticks)
        self.plot_widget.getAxis('left').setTickFont(self.font_ticks)
        self.legend.updateFont(size)

    def mouse_moved(self, evt):
        """Updates the cursor position label when the mouse moves over the plot."""
        pos = evt[0]
        if self.plot_widget.sceneBoundingRect().contains(pos):
            mousePoint = self.plot_widget.plotItem.vb.mapSceneToView(pos)
            self.cursor_pos_label.setText(f"Cursor: x={mousePoint.x():.4f}, y={mousePoint.y():.2f}")

    def updateData(self, data):
        """Public slot to update the plot with new data."""
        self.plot_all_data(data)

    def clear_simulated_data(self):
        """Removes all simulated lines and labels from the plot."""
        while self.simulated_items:
            line, text = self.simulated_items.pop()
            self.plot_widget.removeItem(line)
            self.plot_widget.removeItem(text)
        self.legend.clear()

    def clear_experimental_data(self):
        """Removes the experimental spectrum and peak markers."""
        if self.exp_data_line:
            self.plot_widget.removeItem(self.exp_data_line)
            self.exp_data_line = None
        
        if self.red_triangles:
            self.plot_widget.removeItem(self.red_triangles)
            self.red_triangles = None

    def add_buttons(self, main_layout):
        """Adds the control buttons (Font size, Reset View) to the layout."""
        layout = QHBoxLayout()
        
        font_spin = QSpinBox()
        font_spin.setRange(10, 30)
        font_spin.setValue(self.font_size)
        font_spin.valueChanged.connect(self.update_fonts)
        layout.addWidget(QLabel("Font Size:"))
        layout.addWidget(font_spin)
        
        reset_btn = QPushButton("Reset View")
        reset_btn.clicked.connect(lambda: self.plot_widget.autoRange())
        layout.addWidget(reset_btn)
        
        main_layout.addLayout(layout)