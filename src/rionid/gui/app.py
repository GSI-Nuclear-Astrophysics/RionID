import logging as log
import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QSplitter, QVBoxLayout, QWidget

from .inputs import RionID_GUI
from .plot import CreatePyGUI

log.basicConfig(level=log.DEBUG)


class MainWindow(QWidget):
    """
    The main application window for RionID.

    This class acts as the central container, arranging the input parameters panel
    (left) and the visualization plot (right) using a QSplitter. It handles the
    signal connections between the input logic and the plotting display.

    Attributes
    ----------
    visualization_widget : CreatePyGUI
        The right-hand panel containing the PyQtGraph plot.
    rion_input : RionID_GUI
        The left-hand panel containing input fields and control buttons.
    """

    def __init__(self):
        """
        Initializes the main window, sets geometry, and connects signals.
        """
        super().__init__()
        self.setWindowTitle("RionID")

        screen = QApplication.primaryScreen()
        screen_geom = screen.availableGeometry()
        width = int(screen_geom.width() * 0.8)
        height = int(screen_geom.height() * 0.8)

        # Center the window
        x = (screen_geom.width() - width) // 2
        y = (screen_geom.height() - height) // 2
        self.setGeometry(x, y, width, height)

        # Create a QSplitter to hold both the input and the visualization
        splitter = QSplitter(Qt.Horizontal)

        # Create the visualization before the input panel.
        self.visualization_widget = CreatePyGUI()

        # Create the input widget and connect it to the visualization.
        self.rion_input = RionID_GUI(plot_widget=self.visualization_widget)

        # Add widgets to the splitter
        splitter.addWidget(self.rion_input)
        splitter.addWidget(self.visualization_widget)

        # Set initial size ratios (1 part input, 2 parts plot)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)

        # Create the main layout
        layout = QVBoxLayout()
        layout.addWidget(splitter)
        self.setLayout(layout)

        # --- Signal Connections ---

        # Update plot when a full simulation run finishes
        self.rion_input.visualization_signal.connect(self.update_visualization)

    def update_visualization(self, data):
        """
        Updates the visualization widget with new experimental and simulated data.

        This slot is triggered when a full simulation run is completed. It replaces
        all existing data on the plot (experimental and simulated).

        Parameters
        ----------
        data : ImportData
            The data object containing experimental spectrum arrays and
            simulated ion frequency dictionaries.
        """
        self.visualization_widget.updateData(data)


def main():
    """
    The main entry point for the RionID GUI application.

    Initializes the QApplication, creates the MainWindow, and starts the
    Qt event loop.
    """
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
