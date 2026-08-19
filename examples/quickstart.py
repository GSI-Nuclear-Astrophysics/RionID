"""Minimal, runnable RionID example using synthetic (non-experimental)
data -- no real spectrum or candidate file needed.

Run with: python examples/quickstart.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tests.fixtures.synthetic_spectrum import make_synthetic_spectrum
from rionid.core import ImportData
from rionid.masses import get_ame_data, ionic_moq_u


def main():
    spectrum_path = "quickstart_synthetic.npz"
    make_synthetic_spectrum(spectrum_path)
    print(f"Wrote synthetic spectrum to {spectrum_path}")

    # 72Ge32+ (fully-stripped germanium-72) is the manuscript's own E143
    # worked example (RionID-EPJA/main.tex). A single candidate equal to
    # the reference ion isolates this demo to the frequency model and
    # correction step, without needing a LISE++ candidate file -- see
    # docs/REPRODUCIBILITY.md for a full CLI example with -psim.
    ref_ion = "72Ge32+"
    ame_row = get_ame_data().lookup("Ge", 72)
    ref_moq = ionic_moq_u(ame_row, 32)

    model = ImportData(ref_ion, alphap=0.189, filename=spectrum_path,
                        reload_data=True, circumference=108.36)
    model.moq = {ref_ion: ref_moq}
    model.ref_ion = ref_ion
    model.ref_charge = 32
    model._calculate_srrf(fref=1.93e6)

    print(f"Reference ion {ref_ion}: m/q = {ref_moq} u")
    print(f"Reference frequency: {model.ref_frequency} Hz")
    print(f"srrf (relative revolution frequency): {model.srrf}")
    print("For a full simulation with a real candidate list and CLI/GUI "
          "usage, see docs/REPRODUCIBILITY.md.")


if __name__ == "__main__":
    main()
