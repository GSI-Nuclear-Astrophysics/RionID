import sys
import numpy as np
import re
import pyqtgraph as pg
from PyQt5.QtWidgets import (QMainWindow, QApplication, QVBoxLayout, QWidget, 
                             QPushButton, QHBoxLayout, QLabel, QDesktopWidget, QSpinBox)
from PyQt5.QtGui import QFont, QColor, QBrush
from PyQt5.QtCore import QLoggingCategory, Qt, pyqtSignal

 
class CustomLegendItem(pg.LegendItem):
    """Custom Legend with dynamic font sizing."""
    def __init__(self, font_size, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.font = QFont("Arial", font_size)
        self.brush = pg.mkBrush(255, 255, 255, 200) 
        self.pen = pg.mkPen('k', width=0.5)

    def addItem(self, item, name):
        label = pg.LabelItem(text=name, justify='left')
        label.setFont(self.font)
        super().addItem(item, name)
    
    def updateFont(self, font_size):
        self.font.setPointSize(font_size)
    
    def paint(self, p, *args):
        p.setPen(self.pen)
        p.setBrush(self.brush)
        p.drawRect(self.boundingRect())
        super().paint(p, *args)

class CreatePyGUI(QMainWindow):
    """
    The main visualization widget for RionID.
    """
    plotClicked = pyqtSignal()

    def __init__(self, exp_data=None, sim_data=None):
        super().__init__()
        self.saved_x_range = None  
        self.simulated_items = []
        self.red_triangles = None
        self.exp_data_curve = None
        self.font_size = 14 

        pg.setConfigOptions(antialias=True)  
        pg.setConfigOption('background', 'k') 
        pg.setConfigOption('foreground', 'w')   
        
        self.x_exp = np.array([])
        self.z_exp = np.array([])
        
        self.setup_ui()
        self.plot_widget.scene().sigMouseClicked.connect(self.on_click)

    def on_click(self, event):
        self.plotClicked.emit()

    def setup_ui(self):
        self.setWindowTitle('Schottky Signals Identifier')
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        main_layout = QVBoxLayout(self.main_widget)
        
        QLoggingCategory.setFilterRules('*.warning=false\n*.critical=false')
        
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.showGrid(x=True, y=True, alpha=0.25)
        self.plot_widget.plotItem.ctrl.logYCheck.setChecked(True)
        self.plot_widget.setClipToView(True) 
        
        # Style Axes
        axis_pen = pg.mkPen(color='w', width=1.5)
        self.plot_widget.getAxis('bottom').setPen(axis_pen)
        self.plot_widget.getAxis('left').setPen(axis_pen)
        self.plot_widget.getAxis('bottom').setTextPen('w')
        self.plot_widget.getAxis('left').setTextPen('w')
        
        self.legend = CustomLegendItem(self.font_size, offset=(-10, 10))
        self.legend.brush = pg.mkBrush(0, 0, 0, 150) # Semi-transparent black box
        self.legend.pen = pg.mkPen('w', width=0.5)   # White border
        self.legend.setParentItem(self.plot_widget.graphicsItem())
        
        main_layout.addWidget(self.plot_widget)
        
        self.cursor_pos_label = QLabel(self)
        self.cursor_pos_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.cursor_pos_label)
        self.proxy = pg.SignalProxy(self.plot_widget.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        
        self.add_buttons(main_layout)
        self.update_fonts(self.font_size)

    def _sanitize_positive(self, data, floor=1e-9):
        """
        Aggressively sanitizes data for Log plotting.
        Removes NaNs, Infs, and values <= 0.
        """
        data = np.asanyarray(data, dtype=float)
        # Replace NaNs and Infs with floor
        data[~np.isfinite(data)] = floor
        return np.maximum(data, floor)

    def plot_all_data(self, data):
        # Disable auto-range to prevent ViewBox from calculating bounds on partial data
        self.plot_widget.disableAutoRange()
        try:
            self.clear_experimental_data()
            self.clear_simulated_data()
            self.plot_experimental_data(data)
            self.plot_simulated_data(data)
            
            # Restore view or auto-range if first load
            if self.saved_x_range:
                self.plot_widget.setXRange(*self.saved_x_range, padding=0.02)
            else:
                self.plot_widget.autoRange()
        finally:
            self.plot_widget.enableAutoRange()

    def plot_experimental_data(self, data):
        if data.experimental_data is None: return
        self.exp_data = data.experimental_data
        
        # Extract and Sanitize
        self.x_exp = self.exp_data[0] * 1e-6 # Hz -> MHz
        self.z_exp = self._sanitize_positive(self.exp_data[1])
        
        if len(self.x_exp) == 0: return

        if self.saved_x_range is None:
            self.saved_x_range = (np.min(self.x_exp), np.max(self.x_exp))

        pen = pg.mkPen(color='w', width=1.0)
        brush = pg.mkBrush(color=(255, 255, 255, 50)) # White with low opacity
        
        # Use fillLevel matching the floor to avoid log(-inf) issues
        self.exp_data_curve = pg.PlotCurveItem(
            self.x_exp, self.z_exp, 
            pen=pen, 
            brush=brush, 
            fillLevel=1e-9 
        )
        self.plot_widget.addItem(self.exp_data_curve)
        self.legend.addItem(self.exp_data_curve, 'Experimental Data')
        
        # Plot Peaks
        if hasattr(data, 'peak_freqs') and len(data.peak_freqs) > 0:
            peak_h = self._sanitize_positive(data.peak_heights)
            
            self.red_triangles = self.plot_widget.plot(
                data.peak_freqs * 1e-6, peak_h,
                pen=None, 
                symbol='t1', 
                symbolBrush='#d62728', 
                symbolPen='k',
                symbolSize=10
            )
            self.legend.addItem(self.red_triangles, 'Detected Peaks')

    def plot_simulated_data(self, data):
        self.simulated_data = data.simulated_data_dict
        refion = data.ref_ion
        highlights = data.highlight_ions or []
        
        color_cycle = ['#1f77b4', '#17becf', '#2ca02c', '#d62728', '#9467bd', 
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
        
        color_ref = '#ff7f0e' # Orange for Reference
        color_match = '#2ca02c' # Green for Matches
        
        for i, (harmonic, sdata) in enumerate(self.simulated_data.items()):
            # Arrays for bulk plotting (Vectorization)
            bulk_freqs = []
            bulk_yields = []
            
            # Generate a unique color for this harmonic
            color = color_cycle[i % len(color_cycle)]
            
            for entry in sdata:
                try:
                    freq = float(entry[0]) * 1e-6
                    raw_yield = float(entry[1])
                except (ValueError, TypeError):
                    continue
                
                if not np.isfinite(freq): continue
                yield_value = max(raw_yield, 1.1e-9)
                label = entry[2]
                
                is_highlight = label in highlights
                is_ref = label == refion
                
                # Determine Style
                if is_highlight:
                    c = color_match
                    width = 2
                    style = Qt.SolidLine
                elif is_ref:
                    c = color_ref
                    width = 2
                    style = Qt.DashLine
                else:
                    c = color
                    # No need for width/style here, handled by bulk curve
                
                # Create Text Label
                # Note: Creating thousands of TextItems is still slow. 
                # Use the '-n' (nions) argument to limit this if it's still laggy.
                match = re.match(r'(\d+)([A-Za-z]+)(\d+)\+', label)
                if match:
                    mass, elem, charge = match.groups()
                    new_label = self.to_superscript(mass) + elem + self.to_superscript(charge) + '⁺'
                else: 
                    new_label = label
                
                text_item = pg.TextItem(text=new_label, color=c, anchor=(0.5, 1))
                text_item.setFont(QFont("Arial", self.font_size))
                text_item.setPos(freq, yield_value * 1.05)
                self.plot_widget.addItem(text_item)

                if is_highlight or is_ref:
                    # Plot SPECIAL lines individually (so they draw on top with specific styles)
                    line = self.plot_widget.plot(
                        [freq, freq], [1e-9, yield_value], 
                        pen=pg.mkPen(color=c, width=width, style=style)
                    )
                    self.simulated_items.append((line, text_item))
                else:
                    # Add STANDARD lines to bulk arrays for optimization
                    bulk_freqs.append(freq)
                    bulk_yields.append(yield_value)
                    # Track text item (line is None because it's part of the bulk curve)
                    self.simulated_items.append((None, text_item))

            # --- BULK PLOT ---
            # Draw all standard lines for this harmonic in ONE go
            if bulk_freqs:
                # Interleave arrays for connect='pairs'
                # x: [f1, f1, f2, f2, ...]
                # y: [min, y1, min, y2, ...]
                x_conn = np.repeat(bulk_freqs, 2)
                y_conn = np.empty(len(bulk_yields) * 2)
                y_conn[0::2] = 1e-9
                y_conn[1::2] = bulk_yields
                
                bulk_pen = pg.mkPen(color=color, width=2, style=Qt.DotLine)
                
                # connect='pairs' tells PyQtGraph to draw disjoint lines: (p0->p1), (p2->p3), etc.
                bulk_curve = pg.PlotCurveItem(x_conn, y_conn, connect='pairs', pen=bulk_pen)
                self.plot_widget.addItem(bulk_curve)
                
                # Track the bulk curve so we can clear it later
                self.simulated_items.append((bulk_curve, None))
                
                # Add to legend
                self.legend.addItem(bulk_curve, f'Harmonic {harmonic}')

    def to_superscript(self, s):
        supers = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}
        return ''.join(supers.get(c, c) for c in s)

    def update_fonts(self, size):
        self.font_size = size
        self.font_ticks = QFont("Arial", size)
        self.plot_widget.getAxis('bottom').setTickFont(self.font_ticks)
        self.plot_widget.getAxis('left').setTickFont(self.font_ticks)
        
        label_style = {'color': '#000', 'font-size': f'{size+2}pt'}
        self.plot_widget.setLabel('bottom', 'Frequency (MHz)', **label_style)
        self.plot_widget.setLabel('left', 'Amplitude (a.u.)', **label_style)
        
        self.legend.updateFont(size)

    def mouse_moved(self, evt):
        pos = evt[0]
        if self.plot_widget.sceneBoundingRect().contains(pos):
            mousePoint = self.plot_widget.plotItem.vb.mapSceneToView(pos)
            self.cursor_pos_label.setText(f"Cursor: {mousePoint.x():.4f} MHz")

    def updateData(self, data):
        self.plot_all_data(data)

    def clear_simulated_data(self):
        while self.simulated_items:
            line, text = self.simulated_items.pop()
            if line: 
                self.plot_widget.removeItem(line)
            if text: 
                self.plot_widget.removeItem(text)
        self.legend.clear()

    def clear_experimental_data(self):
        if self.exp_data_curve:
            self.plot_widget.removeItem(self.exp_data_curve)
            self.exp_data_curve = None
        
        if self.red_triangles:
            self.plot_widget.removeItem(self.red_triangles)
            self.red_triangles = None

    def reset_view(self):
        if self.saved_x_range:
            self.plot_widget.setXRange(*self.saved_x_range, padding=0.02)
        
        if len(self.z_exp) > 0:
            min_y = np.min(self.z_exp)
            max_y = np.max(self.z_exp)
            if min_y <= 0: min_y = 1e-9
            self.plot_widget.setYRange(min_y, max_y * 2, padding=0.05)

    def add_buttons(self, main_layout):
        layout = QHBoxLayout()
        
        font_spin = QSpinBox()
        font_spin.setRange(8, 30)
        font_spin.setValue(self.font_size)
        font_spin.valueChanged.connect(self.update_fonts)
        
        lbl = QLabel("Font Size:")
        lbl.setFont(QFont("Arial", 12))
        
        layout.addWidget(lbl)
        layout.addWidget(font_spin)
        
        reset_btn = QPushButton("Reset View")
        reset_btn.setFont(QFont("Arial", 12))
        reset_btn.clicked.connect(self.reset_view)
        layout.addWidget(reset_btn)
        
        main_layout.addLayout(layout)