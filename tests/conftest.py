"""Shared pytest fixtures. QT_QPA_PLATFORM must be set to 'offscreen'
BEFORE any PyQt5 import happens anywhere in the test session -- this
module is pytest's first import, so it goes here.
"""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Preload PyQt5 to prevent pytestqt from trying to detect PySide6
try:
    from PyQt5.QtCore import QCoreApplication  # noqa: F401
except Exception:
    pass

import pytest  # noqa: E402


@pytest.fixture(scope="session")
def qapp():
    """A session-wide QApplication. Runs headlessly via the offscreen
    platform plugin -- confirmed working without Xvfb during the Phase 0
    audit (see docs/PERFORMANCE_BASELINE.md)."""
    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


@pytest.fixture
def synthetic_spectrum_path(tmp_path):
    from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum

    path = tmp_path / "synthetic_spectrum.npz"
    make_synthetic_spectrum(str(path))
    return str(path)
