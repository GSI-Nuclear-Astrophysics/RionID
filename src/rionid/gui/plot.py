import re

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import QLoggingCategory, Qt, QTimer
from PyQt5.QtGui import QClipboard, QFont, QPainterPath
from PyQt5.QtWidgets import (
    QApplication,
    QGraphicsPathItem,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMenu,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from rionid.candidates import matches_selector


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
        self.max_candidate_labels = 300
        self._candidate_data = None
        self._display_options = None
        self._rendering_candidates = False
        self._candidate_timer = QTimer(self)
        self._candidate_timer.setSingleShot(True)
        self._candidate_timer.timeout.connect(self._redraw_candidates)

        pg.setConfigOptions(antialias=True)
        pg.setConfigOption("background", "k")
        pg.setConfigOption("foreground", "w")

        self.x_exp = np.array([])
        self.z_exp = np.array([])

        self.setup_ui()
        self.plot_widget.scene().sigMouseClicked.connect(self._plot_clicked)
        self.plot_widget.sigXRangeChanged.connect(self._schedule_candidate_redraw)

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
        self.candidate_info_label = QLabel("Candidate markers appear after Run", self)
        self.candidate_info_label.setWordWrap(True)
        main_layout.addWidget(self.candidate_info_label)
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
        self._candidate_timer.stop()
        self._display_options = None
        self._rendering_candidates = True
        # Disable auto-range to prevent ViewBox from calculating bounds on partial data
        self.plot_widget.disableAutoRange()
        try:
            self.clear_experimental_data()
            self.clear_simulated_data()
            self.plot_experimental_data(data)
            if self.saved_x_range:
                self.plot_widget.setXRange(*self.saved_x_range, padding=0.02)
            self.plot_simulated_data(data)

            # Restore view or auto-range if first load
            if self.saved_x_range:
                self.plot_widget.setXRange(*self.saved_x_range, padding=0.02)
            else:
                self.plot_widget.autoRange()
        finally:
            self.plot_widget.enableAutoRange()
            self._rendering_candidates = False

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

    def set_display_options(self, display_filter="", highlight_ions="", show_uncertainties=False):
        """Redraw cached candidates without reading inputs or recomputing masses."""
        self._display_options = dict(
            display_filter=display_filter,
            highlight_ions=highlight_ions,
            show_uncertainties=show_uncertainties,
        )
        if self._candidate_data is not None:
            self._redraw_candidates()

    def _schedule_candidate_redraw(self, *_):
        if self._candidate_data is not None and not self._rendering_candidates:
            self._candidate_timer.start(90)

    def _redraw_candidates(self):
        if self._candidate_data is None or self._rendering_candidates:
            return
        self._rendering_candidates = True
        try:
            self.clear_simulated_data()
            if self.exp_data_curve is not None:
                self.legend.addItem(self.exp_data_curve, "Experimental Data")
            self.plot_simulated_data(self._candidate_data)
        finally:
            self._rendering_candidates = False

    def plot_simulated_data(self, data):
        self._candidate_data = data
        options = self._display_options or {}
        selection = options.get("display_filter", getattr(data, "display_filter", ""))
        highlights = options.get("highlight_ions", data.highlight_ions or [])
        show_errors = options.get("show_uncertainties", getattr(data, "show_uncertainties", False))
        source = (
            getattr(data, "all_simulated_data_dict", data.simulated_data_dict)
            if selection
            else data.simulated_data_dict
        )
        self.simulated_data = source
        metadata = getattr(data, "candidate_metadata", {})
        errors = getattr(data, "frequency_uncertainties", {})
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
        label_font = QFont("Arial", self.font_size)
        total = sum(len(rows) for rows in source.values())
        visible_range = self.plot_widget.viewRange()[0]
        if self.saved_x_range is not None and (visible_range[1] - visible_range[0]) < 1e-9:
            visible_range = self.saved_x_range
        margin = max((visible_range[1] - visible_range[0]) * 0.03, 1e-5)
        small_overlay = total <= self.max_candidate_labels
        plotted_count = 0
        labels_left = self.max_candidate_labels
        labels_per_harmonic = max(1, self.max_candidate_labels // max(len(source), 1))
        self.rendered_candidate_names = set()
        self.rendered_label_count = 0
        self.rendered_uncertainty_count = 0
        self.uncertainty_items = []

        for i, (harmonic, rows) in enumerate(source.items()):
            if not len(rows):
                continue
            frequency = np.asarray(rows[:, 0], dtype=float) * 1e-6
            heights = np.asarray(rows[:, 1], dtype=float)
            names = np.asarray(rows[:, 2], dtype=str)
            keep = np.isfinite(frequency)
            if not small_overlay:
                keep &= (frequency >= visible_range[0] - margin) & (
                    frequency <= visible_range[1] + margin
                )
            if selection:
                keep &= np.fromiter(
                    (matches_selector(n, selection) for n in names), bool, len(names)
                )
            indices = np.flatnonzero(keep)
            if not len(indices):
                continue
            color = color_cycle[i % len(color_cycle)]
            positive = heights[np.isfinite(heights) & (heights > 1.1e-9)]
            reference_height = float(np.max(positive)) if len(positive) else 1.0
            unknown = np.fromiter(
                (not metadata.get(n, {}).get("relative_yield_known", True) for n in names),
                bool,
                len(names),
            )
            heights = np.where(np.isfinite(heights), np.maximum(heights, 1.1e-9), 1.1e-9)
            heights[(names == data.ref_ion) & (heights <= 1.1e-9)] = reference_height
            heights[unknown] = max(reference_height * 0.8, 1.1e-9)
            highlighted = np.fromiter(
                (matches_selector(n, highlights) for n in names), bool, len(names)
            )
            reference = names == data.ref_ion
            priority = highlighted | reference | unknown
            ordered = sorted(indices, key=lambda j: (not priority[j], -heights[j], frequency[j]))
            label_indices = set(ordered[: min(labels_left, labels_per_harmonic)])
            labels_left -= len(label_indices)
            plotted_count += len(indices)
            self.rendered_candidate_names.update(names[indices])
            groups = {}
            harmonic_errors = errors.get(str(harmonic), {})
            band_path = QPainterPath()

            for j in indices:
                name, freq, height = str(names[j]), frequency[j], heights[j]
                if highlighted[j]:
                    c, style, width = "#55d96b", Qt.SolidLine, 2
                elif reference[j]:
                    c, style, width = "#ff7f0e", Qt.DashLine, 2
                elif unknown[j]:
                    c, style, width = color, Qt.DashDotLine, 1.8
                else:
                    c, style, width = color, Qt.DotLine, 1.5
                sigma = harmonic_errors.get(name)
                sigma = float(sigma) if sigma is not None else np.nan
                if show_errors and np.isfinite(sigma) and sigma > 0:
                    band_path.addRect(
                        freq - sigma * 1e-6, 1e-9, 2 * sigma * 1e-6, max(height - 1e-9, 1e-9)
                    )
                    self.rendered_uncertainty_count += 1
                text_item = None
                if j in label_indices:
                    match = re.fullmatch(r"(\d+)([A-Za-z]+)(\d+)\+", name)
                    if match:
                        mass, elem, charge = match.groups()
                        label = self.to_superscript(mass) + elem + self.to_superscript(charge) + "⁺"
                    else:
                        label = name
                    text_item = pg.TextItem(text=label, color=c, anchor=(0.5, 1))
                    text_item.setFont(label_font)
                    text_item.setPos(freq, height * 1.05)
                    info = metadata.get(name, {})
                    origin = ", ".join(info.get("sources", [info.get("source", "lise")]))
                    tooltip = f"{name} · harmonic {harmonic}\n{freq * 1e6:.6f} Hz\nSource: {origin}"
                    if unknown[j]:
                        tooltip += "\nNo LISE++ yield; marker height is a display guide."
                    if np.isfinite(sigma):
                        tooltip += f"\nMass-input σf: {sigma:.6g} Hz"
                    if info.get("mass_extrapolated"):
                        tooltip += "\nAME2020 extrapolated mass (#)."
                    text_item.setToolTip(tooltip)
                    self.plot_widget.addItem(text_item)
                    self.rendered_label_count += 1
                if (highlighted[j] or reference[j]) and text_item is not None:
                    line = pg.PlotCurveItem(
                        [freq, freq], [1e-9, height], pen=pg.mkPen(c, width=width, style=style)
                    )
                    self.plot_widget.addItem(line)
                    self.simulated_items.append((line, text_item))
                else:
                    groups.setdefault((c, style, width), []).append((freq, height))
                    if text_item is not None:
                        self.simulated_items.append((None, text_item))
            for (c, style, width), values in groups.items():
                values = np.asarray(values)
                x = np.repeat(values[:, 0], 2)
                y = np.empty(len(x))
                y[0::2], y[1::2] = 1e-9, values[:, 1]
                curve = pg.PlotCurveItem(
                    x, y, connect="pairs", pen=pg.mkPen(c, width=width, style=style)
                )
                self.plot_widget.addItem(curve)
                self.simulated_items.append((curve, None))
            key = pg.PlotCurveItem(pen=pg.mkPen(color, width=2))
            self.legend.addItem(key, f"Harmonic {harmonic}")
            if not band_path.isEmpty():
                band = QGraphicsPathItem(band_path)
                band.setPen(pg.mkPen(None))
                shade = pg.mkColor(color)
                shade.setAlpha(38)
                band.setBrush(pg.mkBrush(shade))
                band.setZValue(-5)
                band.setToolTip(
                    "±1σ propagated AME mass inputs; excludes calibration-model "
                    "and electronic-binding uncertainties."
                )
                self.plot_widget.addItem(band)
                self.simulated_items.append((band, None))
                self.uncertainty_items.append(band)
        self.rendered_candidate_count = plotted_count
        self.candidate_info_label.setText(
            f"{plotted_count:,} candidate markers in view · {self.rendered_label_count:,} labels"
            + (
                " · dash-dot markers: no LISE++ yield"
                if any(not item.get("relative_yield_known", True) for item in metadata.values())
                else ""
            )
        )

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
