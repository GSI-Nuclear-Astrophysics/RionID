"""Golden regression test for the canonical polynomial correction
(ImportData._calculate_srrf + ImportData._simulated_data), cross-checked
against RionID-EPJA/main.tex's own worked numeric example (Eqs. 9-21: the
harmonic-214 -> revolution-space -> harmonic-127 coefficient transform).

See docs/PUBLICATION_TRACEABILITY.md for the manuscript mapping. This test
must keep passing, UNMODIFIED, through every later task in this plan --
none of them may change _calculate_srrf's arithmetic, coefficient
ordering, or units.
"""
from types import SimpleNamespace

import pytest

from rionid.core import ImportData


def test_srrf_correction_matches_manuscript_harmonic_transform():
    # Manuscript's harmonic-214 coefficients (RionID-EPJA/main.tex:372-377)
    a214, b214, c214 = 1.19262945e-9, -9.85167644e-1, 2.03447706e8
    h214 = 214

    # Eq. (14)/(coefficient_h): revolution-space (A, B, C) from a
    # harmonic-h representation.
    A = h214 * a214
    B = b214
    C = c214 / h214
    assert A == pytest.approx(2.55222702e-07, rel=1e-8)
    assert B == pytest.approx(-9.85167644e-01, rel=1e-8)
    assert C == pytest.approx(9.50690215e05, rel=1e-8)

    # Eq. (coefficient_transfer): transform to harmonic 127.
    h127 = 127
    a127 = A / h127
    b127 = B
    c127 = C * h127
    assert a127 == pytest.approx(2.00962758e-09, rel=1e-8)
    assert c127 == pytest.approx(1.20737657e08, rel=1e-8)

    # Now exercise the REAL code path: ImportData._calculate_srrf
    # (Eq. rev_correction/corrected_rev) then ImportData._simulated_data
    # (Eq. harmonic, F=h*f), for a single candidate equal to the reference
    # ion. This isolates srrf to exactly f^(0)/f_ref BEFORE correction
    # (moq[candidate] - moq[ref] == 0), so f^(0) == f_ref exactly and the
    # test can pick f_ref freely without needing a real spectrum, candidate
    # list, or AME table load. See docs/superpowers/specs/
    # 2026-08-19-wave2a-depid-masses-speed-design.md for the derivation.
    model = ImportData.__new__(ImportData)
    model.moq = {"REF+1": 1.0}
    model.ref_ion = "REF+1"
    model.ref_charge = 1
    model.alphap = 0.189
    model.ring = SimpleNamespace(circumference=108.36)
    model.particles_to_simulate = []

    f0_test = 1.0e6  # Hz, arbitrary test revolution frequency
    model._calculate_srrf(fref=f0_test, correct=[A, B, C])

    expected_f_corr = f0_test + (A * f0_test**2 + B * f0_test + C)  # Eq. (corrected_rev)
    actual_f_corr = model.srrf[0] * model.ref_frequency
    assert actual_f_corr == pytest.approx(expected_f_corr, rel=1e-12)

    model._simulated_data(harmonics=[127.0], mode="Frequency")
    actual_F127 = float(model.simulated_data_dict["127.0"][0][0])

    # Direct harmonic-127 polynomial (Eq. h127_corrected), an independent
    # code path from the one above -- this is the actual scaling-invariant
    # check (Eq. scaling_invariant).
    F127_0 = f0_test * h127
    expected_F127 = F127_0 + a127 * F127_0**2 + b127 * F127_0 + c127
    # Matches to 0.000e+00 Hz per the Phase 0 audit; abs tolerance here is
    # generous float64 headroom, not a sign of expected disagreement.
    assert actual_F127 == pytest.approx(expected_F127, abs=1e-6)
