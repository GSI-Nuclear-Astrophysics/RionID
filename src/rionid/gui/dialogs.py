from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QHBoxLayout, 
                             QPushButton, QComboBox)

class KeySelectionDialog(QDialog):
    """
    A modal dialog for selecting data keys from a NumPy (.npz) file.

    When loading generic .npz files, the internal array names are not always
    standard (e.g., they might be 'arr_0', 'arr_1' or named keys like 'frequency').
    This dialog prompts the user to map the available keys to the required
    physics data structures (Frequency and Amplitude).

    Parameters
    ----------
    parent : QWidget, optional
        The parent widget for this dialog.
    keys : list of str, optional
        A list of string keys found in the loaded .npz file.

    Attributes
    ----------
    freq_combo : QComboBox
        Dropdown menu to select the key corresponding to the frequency array.
    amp_combo : QComboBox
        Dropdown menu to select the key corresponding to the amplitude array.
    """

    def __init__(self, parent=None, keys=None):
        """
        Initializes the dialog, sets up the layout, and populates the selection menus.
        """
        super().__init__(parent)
        self.setWindowTitle("Select NPZ Keys")
        
        # Ensure keys is a list
        keys = keys if keys is not None else []

        layout = QVBoxLayout()
        form = QFormLayout()
        
        # Frequency Selection
        self.freq_combo = QComboBox()
        self.freq_combo.addItems(keys)
        
        # Amplitude Selection
        self.amp_combo = QComboBox()
        self.amp_combo.addItems(keys)
        
        # if 'arr_1' exists, set it as default for amplitude
        if len(keys) > 1:
            self.amp_combo.setCurrentIndex(1)

        form.addRow("Frequency Array:", self.freq_combo)
        form.addRow("Amplitude Array:", self.amp_combo)
        
        # Buttons
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        
        layout.addLayout(form)
        layout.addWidget(ok_btn)
        self.setLayout(layout)
        
    def get_params(self):
        """
        Retrieves the user's selection after the dialog is accepted.

        Returns
        -------
        dict
            A dictionary with the following structure:
            {
                'frequency_key': str,
                'amplitude_key': str
            }
        """
        return {
            'frequency_key': self.freq_combo.currentText(),
            'amplitude_key': self.amp_combo.currentText()
        }