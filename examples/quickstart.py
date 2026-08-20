"""Runnable RionID example using synthetic, redistributable inputs.

Run with: python examples/quickstart.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rionid.core import ImportData
from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum


def main():
    data_dir = Path(__file__).resolve().parent / "data"
    data_dir.mkdir(exist_ok=True)
    spectrum_path = data_dir / "quickstart_synthetic.npz"
    make_synthetic_spectrum(spectrum_path)
    print(f"Wrote synthetic spectrum to {spectrum_path}")

    # 72Ge32+ (fully stripped germanium-72) is the worked reference ion.
    ref_ion = "72Ge32+"
    candidates_path = data_dir / "candidates.lpp"

    model = ImportData(
        ref_ion, alphap=0.189, filename=spectrum_path, reload_data=True, circumference=108.36
    )
    model._set_particles_to_simulate_from_file(candidates_path)
    model._calculate_moqs()
    model._calculate_srrf(fref=1.93e6)
    model._simulated_data(harmonics=[127.0], mode="frequency")

    print(f"Loaded {len(model.nuclei_names)} candidates from {candidates_path}")
    print(f"Reference ion {ref_ion}: m/q = {model.moq[ref_ion]} u")
    print(f"Reference frequency: {model.ref_frequency} Hz")
    print(f"Simulated harmonic 127 shape: {model.simulated_data_dict['127.0'].shape}")


if __name__ == "__main__":
    main()
