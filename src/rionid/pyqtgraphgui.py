import sys
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QLabel, QDesktopWidget
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QLoggingCategory, Qt

class CreatePyGUI(QMainWindow):
    '''
    PyView (MVC)
    '''
    def __init__(self):
        super().__init__()
        self.saved_x_range = None  
        self.saved_y_range = None 
        self.simulated_items = []
        self.setup_ui()

    def setup_ui(self):
        self.setWindowTitle('Schottky Signals Identifier')
        width = QDesktopWidget().screenGeometry(-1).width()
        height = QDesktopWidget().screenGeometry(-1).height()
        self.setGeometry(100, 100, width, height)  # Set window size
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        main_layout = QVBoxLayout(self.main_widget)
        #logging annoying messages
        QLoggingCategory.setFilterRules('*.warning=false\n*.critical=false')
        # Create the plot widget
        self.plot_widget = pg.PlotWidget()
        
        # Set logY to true
        self.plot_widget.plotItem.ctrl.logYCheck.setChecked(True)

        main_layout.addWidget(self.plot_widget)
        # Add legend
        self.legend = pg.LegendItem(offset=(-10, 10))  # Adjust offset as needed
        self.legend.setParentItem(self.plot_widget.graphicsItem())
        self.plot_widget.setLabel(
            "left",
            '<span style="color: white; font-size: 20px">Amplitude (arb. units)</span>'
        )
        self.plot_widget.setLabel(
            "bottom",
            '<span style="color: white; font-size: 20px">Frequency (MHz)</span>'
        )
        
        # Customizing tick label font size
        font_ticks = QFont()
        font_ticks.setPixelSize(20)  # Set the desired font size here
        self.plot_widget.getAxis('bottom').setTickFont(font_ticks)
        self.plot_widget.getAxis("bottom").setStyle(tickTextOffset = 15)
        self.plot_widget.getAxis('left').setTickFont(font_ticks)
        self.plot_widget.getAxis("left").setStyle(tickTextOffset = 15)
        
        # Cursor position label
        font = QFont("Times", 12)
        self.cursor_pos_label = QLabel(self)
        self.cursor_pos_label.setFont(font)
        main_layout.addWidget(self.cursor_pos_label)
        self.proxy = pg.SignalProxy(self.plot_widget.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        
        # Add control buttons
        self.add_buttons(main_layout)
    
    def plot_all_data(self, data):
        self.plot_widget.clear()
        self.plot_experimental_data(data)
        self.plot_simulated_data(data)

    def plot_experimental_data(self, data):
        if data.experimental_data is None:  # Check if experimental data is available
            print("No experimental data available, skipping experimental data plotting.")
            return  # Skip plotting experimental data
        self.exp_data = data.experimental_data
        # Plot experimental data
        if hasattr(self, 'exp_data_line'):
            self.plot_widget.removeItem(self.exp_data_line)

        # Set the initial X-range to encompass all experimental data
        self.x_exp, self.z_exp = self.exp_data[0]*1e-6, self.exp_data[1]
        
        if self.saved_x_range is None:
            print(" 1.self.saved_x_range is None: ")
            self.saved_x_range = (min(self.x_exp), max(self.x_exp))
            self.plot_widget.setXRange(*self.saved_x_range, padding=0.05)

            # Save the Y range as well
            min_z = np.min(self.z_exp)
            max_z = max(self.z_exp)
            if min_z <= 0:
                # Handle the logarithmic scale by setting the minimum to a small value if necessary
                min_z = 1e-10  # or some other small positive value
            print("min_z ",min_z ," max_z = ",max_z)

        self.exp_data_line = self.plot_widget.plot(self.x_exp, self.z_exp, pen=pg.mkPen('white', width=3))
        self.legend.addItem(self.exp_data_line, 'Experimental Data')

    def plot_simulated_data(self, data):
        self.simulated_data = data.simulated_data_dict
        refion = data.ref_ion
        for i, (harmonic, sdata) in enumerate(self.simulated_data.items()):
            color = pg.intColor(i, hues=len(self.simulated_data))
            for entry in sdata:
                freq = float(entry[0])*1e-6
                label = entry[2]
                yield_value = float(entry[1])  # Simulated yield (amplitude)
                #print("freq = ",freq," yield_value = ",yield_value)
                # Find the corresponding z_exp for the given freq
                freq_range = 0.005
                if data.experimental_data is None:  # Check if experimental data is available
                    z_value = yield_value
                else:
                    z_value = self.get_z_exp_at_freq(freq, freq_range)
                label_color = None
                # Set label color to yellow if it matches the reference ion
                if label == refion:
                    label_color = 'yellow'  # If matching, use yellow
                else:
                    label_color = color  # Otherwise, use the default color
                
                # Vertical line
                line = self.plot_widget.plot([freq, freq], [1e-10, z_value], pen=pg.mkPen(color=label_color, width=1, style = Qt.DashLine ))
                # Text label at top
                text = pg.TextItem(text=label, color=label_color, anchor=(0.5, 0))
                self.plot_widget.addItem(text)

                logy_checked = self.plot_widget.plotItem.ctrl.logYCheck.isChecked()
                if logy_checked:
                    text.setPos(freq, np.log10(z_value) + 0.2)  # Adjust 0.1 as needed for visibility
                else:
                    text.setPos(freq, z_value * 1.05)
                self.simulated_items.append((line, text))  # Add as a tuple

            self.legend.addItem(line, f'Harmonic = {float(harmonic)} ; Bρ = {data.brho:.6f} [Tm].')
            
    def get_z_exp_at_freq(self, freq, freq_range):
        # Check if self.x_exp and self.z_exp are not empty
        if len(self.x_exp) == 0 or len(self.z_exp) == 0:
            return None

        # Define the frequency range
        lower_bound = freq - freq_range
        upper_bound = freq + freq_range
        
        # Find indices within the specified frequency range
        indices = (self.x_exp >= lower_bound) & (self.x_exp <= upper_bound)
        if not np.any(indices):
            return None
        
        # Return the maximum z_exp value within the specified range
        return np.max(self.z_exp[indices])
                                                        
    def updateData(self, data):
        print("Updating data in visualization GUI...")
        self.clear_experimental_data()
        self.clear_simulated_data()
        self.plot_all_data(data)

    def clear_simulated_data(self):
        print("Clearing simulated data plots...")
        while self.simulated_items:
            line, text = self.simulated_items.pop()
            self.plot_widget.removeItem(line)
            self.plot_widget.removeItem(text)
        self.legend.clear()
        self.simulated_data = None

    def clear_experimental_data(self):
        if hasattr(self, 'exp_data_line'):
            print("Clearing experimental data plot...")
            self.plot_widget.removeItem(self.exp_data_line)
            self.legend.removeItem(self.exp_data_line)
        self.exp_data_line = None
        self.exp_data = None

    def toggle_simulated_data(self):
        for line, text in self.simulated_items:
            line.setVisible(not line.isVisible())
            text.setVisible(not text.isVisible())

    def mouse_moved(self, evt):
        pos = evt[0]  # using signal proxy turns original arguments into a tuple
        if self.plot_widget.sceneBoundingRect().contains(pos):
            mousePoint = self.plot_widget.plotItem.vb.mapSceneToView(pos)
            self.cursor_pos_label.setText(f"Cursor Position: x={mousePoint.x():.8f}, y={mousePoint.y():.2f}")
            # Check the state of the LogY checkbox
            #logy_checked = self.plot_widget.plotItem.ctrl.logYCheck.isChecked()
            #print(f"LogY is {'enabled' if logy_checked else 'disabled'}")
                                    
    def save_selected_data(self):
        selected_range = self.plot_widget.getViewBox().viewRange()[0]
        mask = (self.x_exp >= selected_range[0]) & (self.x_exp <= selected_range[1])
        selected_data = self.z_exp[mask]
        selected_x = self.x_exp[mask]
        filename = 'selected_data.npz'
        np.savez(filename, x=selected_x, z=selected_data)
        print(f"Data saved to {filename}")

    def reset_view(self):
        # Reset the plot to the original X and Y ranges
        self.plot_widget.setXRange(*self.initial_x_range, padding=0.05)
        self.plot_widget.setYRange(min(self.z_exp), max(self.z_exp), padding=0.05)

    def add_buttons(self, main_layout):
        button_layout = QHBoxLayout()

        font = QFont("Times", 15)
        font.setBold(True)

        toggle_button = QPushButton("Toggle Simulated Data")
        toggle_button.clicked.connect(self.toggle_simulated_data)
        toggle_button.setFont(font)
        button_layout.addWidget(toggle_button)

        save_button = QPushButton("Save Selected Data")
        save_button.clicked.connect(self.save_selected_data)
        save_button.setFont(font)
        button_layout.addWidget(save_button)

        zoom_in_button = QPushButton("Zoom In")
        zoom_in_button.clicked.connect(lambda: self.plot_widget.getViewBox().scaleBy((0.5, 0.5)))
        zoom_in_button.setFont(font)
        button_layout.addWidget(zoom_in_button)

        zoom_out_button = QPushButton("Zoom Out")
        zoom_out_button.clicked.connect(lambda: self.plot_widget.getViewBox().scaleBy((2, 2)))
        zoom_out_button.setFont(font)
        button_layout.addWidget(zoom_out_button)

        # Button to reset the plot view
        reset_view_button = QPushButton("Reset View")
        reset_view_button.setFont(font)
        reset_view_button.clicked.connect(self.reset_view)
        button_layout.addWidget(reset_view_button)

        main_layout.addLayout(button_layout)


# Example Usage:
if __name__ == '__main__':
    app = QApplication(sys.argv)

    # Example data
    experimental_data = np.array([[2.35000019e+08, 9.04612897e-02],
                                  [2.35000057e+08, 9.07298288e-02],
                                  [2.35000095e+08, 9.01448335e-02],
                                  [2.54999905e+08, 9.01264557e-02],
                                  [2.54999943e+08, 9.01772547e-02],
                                  [2.54999981e+08, 9.03425368e-02]])

    simulated_data = {
        '1.0': np.array([['241127381.22165576', '0.00054777', '80Kr+35'],
                         ['242703150.0762615', '0.0048654', '79Br+35']])
    }

    sa = CreatePyGUI(experimental_data, simulated_data)
    sa.show()
    sys.exit(app.exec_())
