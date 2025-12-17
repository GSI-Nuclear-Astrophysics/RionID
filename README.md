# RionID (Ring-stored ion IDentification)

[![Documentation](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://GSI-Nuclear-Astrophysics.github.io/RionID)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8169341.svg)](https://doi.org/10.5281/zenodo.8169341)
[![PyPI version](https://badge.fury.io/py/rionid.svg)](https://badge.fury.io/py/rionid)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**RionID** is a Python software for the identification of ions stored in storage rings. It simulates revolution frequencies based on magnetic rigidity or frequency settings and matches them against experimental Schottky spectra.

<div align="center">
  <img src="https://github.com/GSI-Nuclear-Astrophysics/rionid/raw/master/docs/img/rionid.png?raw=true" width="50%">
</div>

## Features
*   **Pure Python:** No ROOT dependencies required.
*   **Automated Matching:** Includes Quick PID logic to scan $\alpha_p$ and Reference Frequency to find the best match ($\chi^2$ minimization).
*   **Signal Processing:** Built-in baseline subtraction (BrPLS) and peak detection.
*   **Standalone:** Bundles `barion` and `lisereader` (GPL-3.0) for easy installation without complex dependency management.

## Installation

### Option 1: From PyPI (Recommended)
RionID is available on the Python Package Index. This is the easiest way to install it along with all dependencies.

```bash
pip install rionid
```

### Option 2: From Source
If you want the development version:

```bash
git clone https://github.com/GSI-Nuclear-Astrophysics/rionid.git
cd rionid
pip install .
```

## Usage

### Graphical User Interface (GUI)
Once installed, you can launch the GUI simply by typing:

```bash
rionid
```

### Command Line Interface (CLI)
You can also run simulations directly from the terminal:

```bash
rionid datafile.npz -f 11.2452 -r 209Bi+83 -psim fragments.lpp -b 5.5
```

## Arguments

*   `datafile`: Input spectrum file (.npz, .csv, .txt).
*   `-r`, `--refion`: Reference ion (e.g., `72Ge+35`).
*   `-ap`, `--alphap`: Momentum compaction factor.
*   `-psim`: LISE++ output file for fragment yields.
*   `-hrm`: Harmonics to simulate.
*   `-b`, `--brho`: Magnetic rigidity (Brho) [Tm].
*   `-f`, `--fref`: Revolution frequency [Hz].
*   `--remove_baseline`: Apply baseline subtraction.
*   `--peak_threshold_pct`: Peak detection threshold (0.0 - 1.0).

## Acknowledgements
*   **Dr. RuiJiu Chen** for providing the C++ Time-of-Flight simulation code that inspired the backbone of this software.
*   **Dr. Shahab Sanjari** for guidance on software architecture and Schottky analysis.

## License
This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.
```