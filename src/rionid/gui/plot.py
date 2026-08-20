import re

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import QLoggingCategory, Qt
from PyQt5.QtGui import QClipboard, QFont
from PyQt5.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMenu,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class CustomLegendItem(pg.LegendItem):
    """Custom Legend with dynamic font sizing."""

    def __init__(self, font_size, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.font = QFont("Arial", font_size)
        self.brush = pg.mkBrush(255, 255, 255, 200)
        self.pen = pg.mkPen("k", width=0.5)

    def addItem(self, item, name):
        label = pg.LabelItem(text=name, justify="left")
        label.setFont(self.font)
        super().addItem(item, name)

    def updateFont(self, font_size):
        self.font.setPointSize(font_size)

    def paint(self, p, *args):
        p.setPen(self.pen)
        p.setBrush(self.brush)
        p.drawRect(self.boundingRect())
        super().paint(p, *args)


class FrequencyMarker(pg.InfiniteLine):
    """Movable frequency marker with copy/remove actions."""

    def __init__(
        self,
        marker_id,
        frequency_mhz,
        remove_callback,
        selection_callback,
        position_callback,
        copied_callback,
    ):
        super().__init__(
            pos=frequency_mhz,
            angle=90,
            movable=True,
            pen=pg.mkPen("#ffd700", width=2),
            hoverPen=pg.mkPen("#ffffff", width=3),
        )
        self.marker_id = marker_id
        self.remove_callback = remove_callback
        self.selection_callback = selection_callback
        self.position_callback = position_callback
        self.copied_callback = copied_callback
        self.selected = False
        self.frequency_label = pg.InfLineLabel(
            self, self._label_text(), position=0.95, color="#ffd700"
        )
        self.sigPositionChanged.connect(self._position_changed)

    def _label_text(self):
        return f"M{self.marker_id}: {float(self.value()):.9f} MHz"

    def _position_changed(self):
        self.frequency_label.setText(self._label_text())

        self.position_callback()

    def set_selected(self, selected):
        self.selected = selected
        color = "#00e5ff" if selected else "#ffd700"
        self.setPen(pg.mkPen(color, width=3 if selected else 2))
        self.frequency_label.setColor(color)

    def copy_frequency_hz(self):
        text = f"{float(self.value()) * 1e6:.3f} Hz"
        QApplication.clipboard().setText(text, QClipboard.Clipboard)
        QApplication.processEvents()
        self.copied_callback(text)
        return text

    def mouseClickEvent(self, ev):
        if ev.button() == Qt.LeftButton:
            self.selection_callback(self)
            super().mouseClickEvent(ev)
            return
        if ev.button() != Qt.RightButton:
            super().mouseClickEvent(ev)
            return

        menu = QMenu()
        copy_action = menu.addAction("Copy frequency (Hz)")
        remove_action = menu.addAction("Remove marker")
        selected = menu.exec_(ev.screenPos().toPoint())
        if selected == copy_action:
            self.copy_frequency_hz()
        elif selected == remove_action:
            self.remove_callback(self)
        ev.accept()


class CreatePyGUI(QMainWindow):
    """
    The main visualization widget for RionID.
    """

    def __init__(self, exp_data=None, sim_data=None):
        super().__init__()
        self.saved_x_range = None
        self.simulated_items = []
        self.frequency_markers = []
        self.selected_markers = []
        self.next_marker_id = 1
        self.exp_data_curve = None
        self.font_size = 14

        pg.setConfigOptions(antialias=True)
        pg.setConfigOption("background", "k")
        pg.setConfigOption("foreground", "w")

        self.x_exp = np.array([])
        self.z_exp = np.array([])

        self.setup_ui()
        self.plot_widget.scene().sigMouseClicked.connect(self._plot_clicked)

    def setup_ui(self):
        self.setWindowTitle("Schottky Signals Identifier")
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        main_layout = QVBoxLayout(self.main_widget)

        QLoggingCategory.setFilterRules("*.warning=false\n*.critical=false")

        self.plot_widget = pg.PlotWidget()
        self.plot_widget.showGrid(x=True, y=True, alpha=0.25)
        self.plot_widget.plotItem.ctrl.logYCheck.setChecked(True)
        self.plot_widget.setClipToView(True)

        # Style Axes
        axis_pen = pg.mkPen(color="w", width=1.5)
        self.plot_widget.getAxis("bottom").setPen(axis_pen)
        self.plot_widget.getAxis("left").setPen(axis_pen)
        self.plot_widget.getAxis("bottom").setTextPen("w")
        self.plot_widget.getAxis("left").setTextPen("w")

        self.legend = CustomLegendItem(self.font_size, offset=(-10, 10))
        self.legend.brush = pg.mkBrush(0, 0, 0, 150)  # Semi-transparent black box
        self.legend.pen = pg.mkPen("w", width=0.5)  # White border
        self.legend.setParentItem(self.plot_widget.graphicsItem())

        main_layout.addWidget(self.plot_widget)

        self.cursor_pos_label = QLabel(self)
        self.cursor_pos_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.cursor_pos_label)
        self.marker_info_label = QLabel("Select two markers to measure Δf", self)
        self.marker_info_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.marker_info_label)
        self.proxy = pg.SignalProxy(
            self.plot_widget.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved
        )

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
        if data.experimental_data is None:
            return
        self.exp_data = data.experimental_data

        # Extract and Sanitize
        self.x_exp = self.exp_data[0] * 1e-6  # Hz -> MHz
        self.z_exp = self._sanitize_positive(self.exp_data[1])

        if len(self.x_exp) == 0:
            return

        if self.saved_x_range is None:
            self.saved_x_range = (np.min(self.x_exp), np.max(self.x_exp))

        pen = pg.mkPen(color="w", width=1.0)
        brush = pg.mkBrush(color=(255, 255, 255, 50))  # White with low opacity

        # Use fillLevel matching the floor to avoid log(-inf) issues
        self.exp_data_curve = pg.PlotCurveItem(
            self.x_exp, self.z_exp, pen=pen, brush=brush, fillLevel=1e-9
        )
        self.plot_widget.addItem(self.exp_data_curve)
        self.legend.addItem(self.exp_data_curve, "Experimental Data")

    def plot_simulated_data(self, data):
        # Profiling at N=2000 candidates shows this method dominated not by
        # pg.TextItem construction itself (~0.14ms/item) but by
        # PyQtGraph's addItem/removeItem scene-graph reparenting overhead
        # (itemChange/changeParent/signal connect-disconnect), triggered
        # on every clear_simulated_data()+plot_simulated_data() cycle --
        # roughly 0.7-1.2s of the ~1.4s total at N=2000. The real fix is
        # reusing TextItem objects across redraws instead of destroying
        # and recreating them, which requires diffing the previous vs.
        # new (harmonic, ion) label set and correctly handling ions that
        # change highlight/reference status between redraws (their
        # colour must update, not just position/text). That correctness
        # needs visual QA before implementation.
        self.simulated_data = data.simulated_data_dict
        refion = data.ref_ion
        highlights = data.highlight_ions or []

        color_cycle = [
            "#1f77b4",
            "#17becf",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
        ]

        color_ref = "#ff7f0e"  # Orange for Reference
        color_match = "#2ca02c"  # Green for Matches
        label_font = QFont("Arial", self.font_size)

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

                if not np.isfinite(freq):
                    continue
                yield_value = max(raw_yield, 1.1e-9)
                label = entry[2]

                is_highlight = label in highlights
                is_ref = label == refion

                # Determine Style
                width = style = None  # only read under is_highlight/is_ref, set below
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
                match = re.match(r"(\d+)([A-Za-z]+)(\d+)\+", label)
                if match:
                    mass, elem, charge = match.groups()
                    new_label = self.to_superscript(mass) + elem + self.to_superscript(charge) + "⁺"
                else:
                    new_label = label

                text_item = pg.TextItem(text=new_label, color=c, anchor=(0.5, 1))
                text_item.setFont(label_font)
                text_item.setPos(freq, yield_value * 1.05)
                self.plot_widget.addItem(text_item)

                if is_highlight or is_ref:
                    # Plot SPECIAL lines individually (so they draw on top with specific styles)
                    line = self.plot_widget.plot(
                        [freq, freq],
                        [1e-9, yield_value],
                        pen=pg.mkPen(color=c, width=width, style=style),
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
                bulk_curve = pg.PlotCurveItem(x_conn, y_conn, connect="pairs", pen=bulk_pen)
                self.plot_widget.addItem(bulk_curve)

                # Track the bulk curve so we can clear it later
                self.simulated_items.append((bulk_curve, None))

                # Add to legend
                self.legend.addItem(bulk_curve, f"Harmonic {harmonic}")

    def to_superscript(self, s):
        supers = {
            "0": "⁰",
            "1": "¹",
            "2": "²",
            "3": "³",
            "4": "⁴",
            "5": "⁵",
            "6": "⁶",
            "7": "⁷",
            "8": "⁸",
            "9": "⁹",
        }
        return "".join(supers.get(c, c) for c in s)

    def update_fonts(self, size):
        self.font_size = size
        self.font_ticks = QFont("Arial", size)
        self.plot_widget.getAxis("bottom").setTickFont(self.font_ticks)
        self.plot_widget.getAxis("left").setTickFont(self.font_ticks)

        label_style = {"color": "#000", "font-size": f"{size + 2}pt"}
        self.plot_widget.setLabel("bottom", "Frequency (MHz)", **label_style)
        self.plot_widget.setLabel("left", "Amplitude (a.u.)", **label_style)

        self.legend.updateFont(size)

    def mouse_moved(self, evt):
        pos = evt[0]
        if self.plot_widget.sceneBoundingRect().contains(pos):
            mousePoint = self.plot_widget.plotItem.vb.mapSceneToView(pos)
            self.cursor_pos_label.setText(self._format_cursor_position(mousePoint.x()))

    @staticmethod
    def _format_cursor_position(frequency_mhz):
        """Format MHz with enough decimal places to resolve millihertz."""
        return f"Cursor: {frequency_mhz:.9f} MHz"

    def _plot_clicked(self, event):
        if event.button() != Qt.LeftButton or not event.double():
            return
        if not self.plot_widget.sceneBoundingRect().contains(event.scenePos()):
            return
        point = self.plot_widget.plotItem.vb.mapSceneToView(event.scenePos())
        self.add_frequency_marker(point.x())
        event.accept()

    def add_frequency_marker(self, frequency_mhz):
        marker = FrequencyMarker(
            self.next_marker_id,
            frequency_mhz,
            self.remove_frequency_marker,
            self.toggle_marker_selection,
            self.update_marker_measurement,
            self.show_copy_confirmation,
        )
        self.next_marker_id += 1
        self.frequency_markers.append(marker)
        self.plot_widget.addItem(marker)
        return marker

    def remove_frequency_marker(self, marker):
        if marker in self.frequency_markers:
            if marker in self.selected_markers:
                self.selected_markers.remove(marker)
            self.frequency_markers.remove(marker)
            self.plot_widget.removeItem(marker)
            self.update_marker_measurement()

    def clear_frequency_markers(self):
        for marker in list(self.frequency_markers):
            self.remove_frequency_marker(marker)

    def toggle_marker_selection(self, marker):
        if marker in self.selected_markers:
            self.selected_markers.remove(marker)
            marker.set_selected(False)
        else:
            if len(self.selected_markers) == 2:
                oldest = self.selected_markers.pop(0)
                oldest.set_selected(False)
            self.selected_markers.append(marker)
            marker.set_selected(True)
        self.update_marker_measurement()

    def update_marker_measurement(self):
        if len(self.selected_markers) != 2:
            self.marker_info_label.setText("Select two markers to measure Δf")
            return
        first, second = self.selected_markers
        delta_hz = abs(float(second.value()) - float(first.value())) * 1e6
        self.marker_info_label.setText(
            f"Δf M{first.marker_id} ↔ M{second.marker_id}: {delta_hz:.3f} Hz"
        )

    def copy_selected_markers(self):
        if not self.selected_markers:
            self.marker_info_label.setText("Select a marker before copying")
            return ""
        text = "\n".join(
            f"{float(marker.value()) * 1e6:.3f} Hz" for marker in self.selected_markers
        )
        QApplication.clipboard().setText(text, QClipboard.Clipboard)
        QApplication.processEvents()
        self.show_copy_confirmation(text)
        return text

    def show_copy_confirmation(self, text):
        self.marker_info_label.setText(f"Copied: {text.replace(chr(10), ', ')}")

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

    def reset_view(self):
        if self.saved_x_range:
            self.plot_widget.setXRange(*self.saved_x_range, padding=0.02)

        if len(self.z_exp) > 0:
            min_y = np.min(self.z_exp)
            max_y = np.max(self.z_exp)
            if min_y <= 0:
                min_y = 1e-9
            self.plot_widget.setYRange(min_y, max_y * 2, padding=0.05)

    def add_buttons(self, main_layout):
        layout = QHBoxLayout()

        self.reset_view_button = QPushButton("Reset View")
        self.reset_view_button.setFont(QFont("Arial", 12))
        self.reset_view_button.setStyleSheet(
            "QPushButton { background-color: #bbdefb; color: black; }"
        )
        self.reset_view_button.clicked.connect(self.reset_view)
        layout.addWidget(self.reset_view_button)

        self.copy_markers_button = QPushButton("Copy Selected")
        self.copy_markers_button.setFont(QFont("Arial", 12))
        self.copy_markers_button.setStyleSheet(
            "QPushButton { background-color: #c8e6c9; color: black; }"
        )
        self.copy_markers_button.clicked.connect(self.copy_selected_markers)
        layout.addWidget(self.copy_markers_button)

        self.clear_markers_button = QPushButton("Clear Markers")
        self.clear_markers_button.setFont(QFont("Arial", 12))
        self.clear_markers_button.setStyleSheet(
            "QPushButton { background-color: #fff9c4; color: black; }"
        )
        self.clear_markers_button.clicked.connect(self.clear_frequency_markers)
        layout.addWidget(self.clear_markers_button)

        main_layout.addLayout(layout)
