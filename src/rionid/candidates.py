"""Candidate selection and cached AME2020 mass-to-charge searches."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .masses import get_ame_data, ionic_mass_u

ELEMENTS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs "
    "Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl "
    "Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh "
    "Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
ATOMIC_NUMBER = {element: z for z, element in enumerate(ELEMENTS)}
ISOMER_LABELS = "mnpqrx"


@dataclass(frozen=True)
class Ion:
    A: int
    element: str
    q: int
    isomer: str = ""

    @property
    def name(self):
        return f"{self.A}{self.element}{self.isomer}{self.q}+"

    @property
    def Z(self):
        return ATOMIC_NUMBER[self.element]


def split_ions(text):
    if not text:
        return []
    return [x for x in re.split(r"[\s,;]+", text.strip()) if x]


@lru_cache(maxsize=8192)
def parse_ion(text: str) -> Ion:
    """Accept 205Au76+ and 205Au+76, requiring a physical positive charge.

    A trailing Nubase isomer label (m, n, p, q, r, x) names an excited state,
    as in 180Hfm76+.
    """
    match = re.fullmatch(r"(\d+)([A-Za-z]{1,3})(?:(\d+)\+|\+(\d+))", text.strip())
    if not match:
        raise ValueError(f"Invalid ion {text!r}; use an explicit charge, for example 205Au76+")
    mass, letters, q1, q2 = match.groups()
    q = int(q1 or q2)
    readings = [(letters.capitalize(), "")]
    if len(letters) > 1 and letters[-1].lower() in ISOMER_LABELS:
        readings.append((letters[:-1].capitalize(), letters[-1].lower()))
    known = [(e, i) for e, i in readings if e in ATOMIC_NUMBER and e != "n"]
    if not known:
        raise ValueError(f"Unknown element in {text!r}")
    for element, isomer in known:
        if int(mass) >= ATOMIC_NUMBER[element] and 1 <= q <= ATOMIC_NUMBER[element]:
            return Ion(int(mass), element, q, isomer)
    raise ValueError(f"Invalid mass number or charge in {text!r}; require 1 <= q <= Z")


@lru_cache(maxsize=1024)
def _selectors(text):
    compiled = []
    for token in split_ions(text):
        token = token.rstrip("*")
        if "+" in token:
            compiled.append(("ion", parse_ion(token).name))
            continue
        match = re.fullmatch(r"(\d+)?([A-Za-z]{1,2})", token)
        if not match:
            raise ValueError(f"Invalid selector {token!r}; use 205Au76+, 205Au or Au")
        mass, element = match.groups()
        element = element.capitalize()
        if element not in ATOMIC_NUMBER:
            raise ValueError(f"Unknown element in selector {token!r}")
        compiled.append(
            ("nuclide" if mass else "element", (int(mass), element) if mass else element)
        )
    return tuple(compiled)


def matches_selector(ion_name, selectors):
    """Match exact ions, an isotope in every charge state, or an element."""
    if not selectors:
        return False
    text = selectors if isinstance(selectors, str) else " ".join(selectors)
    ion = parse_ion(str(ion_name))
    return any(
        (kind == "ion" and value == ion.name)
        or (kind == "nuclide" and value == (ion.A, ion.element))
        or (kind == "element" and value == ion.element)
        for kind, value in _selectors(text)
    )


def validate_selectors(text):
    _selectors(text)


def mass_metadata(row, source):
    value = row[9]
    sigma = float(value) if value is not None and math.isfinite(float(value)) else None
    return {
        "source": source,
        "sources": [source],
        "mass_uncertainty_keV": sigma,
        "mass_extrapolated": str(row[0]).strip() == "sys",
        "relative_yield_known": source == "lise",
    }


def ion_from_ame(name, ame=None):
    ion = parse_ion(name)
    if ion.isomer:
        raise ValueError(f"{ion.name} is an isomer; use the Show isomers option")
    ame = ame or get_ame_data()
    row = ame.lookup(ion.element, ion.A)
    if row is None:
        raise ValueError(f"{ion.A}{ion.element} has no mass in AME2020")
    mass = float(ionic_mass_u(row, ion.q))
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError(f"{ion.name} has no usable AME2020 mass")
    return ion, row, mass


_CATALOGUES = {}


def _catalogue(ame):
    """Build and sort each table once; searches use binary interval lookup."""
    cached = _CATALOGUES.get(id(ame))
    if cached is not None and cached[0] is ame:
        return cached[1:]
    names, moqs, rows, charges = [], [], [], []
    seen = set()
    for row in ame.ame_table:
        z, mass_number, element = int(row[4]), int(row[5]), row[6]
        if (element, mass_number) in seen or not 1 <= z < len(ELEMENTS):
            continue
        seen.add((element, mass_number))
        for q in range(1, z + 1):
            mass = float(ionic_mass_u(row, q))
            if not math.isfinite(mass) or mass <= 0:
                continue
            names.append(f"{mass_number}{element}{q}+")
            moqs.append(mass / q)
            rows.append(row)
            charges.append(q)
    order = np.argsort(moqs)
    moqs = np.asarray(moqs)[order]
    names = np.asarray(names)[order]
    charges = np.asarray(charges)[order]
    rows = [rows[i] for i in order]
    _CATALOGUES[id(ame)] = (ame, moqs, names, charges, rows)
    return moqs, names, charges, rows


def search_moq(center_ion, tolerance_percent=1.0, charge_min=1, charge_max=None, ame=None):
    """Return every AME2020 ionic combination within an inclusive m/q interval.

    Search includes AME extrapolations and labels them in the returned metadata.
    Entries are candidate ground states, without a production-yield assumption.
    """
    if not math.isfinite(tolerance_percent) or not 0 <= tolerance_percent < 100:
        raise ValueError("Mass-to-charge tolerance must be finite and in [0,100) percent")
    if not math.isfinite(charge_min) or charge_min < 1 or int(charge_min) != charge_min:
        raise ValueError("Minimum charge must be a positive integer")
    if charge_max is not None and (
        not math.isfinite(charge_max) or int(charge_max) != charge_max or charge_max < charge_min
    ):
        raise ValueError("Maximum charge must be an integer >= minimum charge")
    ame = ame or get_ame_data()
    center, _, mass = ion_from_ame(center_ion, ame)
    target = mass / center.q
    moqs, names, charges, rows = _catalogue(ame)
    low, high = target * (1 - tolerance_percent / 100), target * (1 + tolerance_percent / 100)
    first = np.searchsorted(moqs, low, side="left")
    last = np.searchsorted(moqs, high, side="right")
    chosen = np.arange(first, last)
    mask = charges[chosen] >= charge_min
    if charge_max is not None:
        mask &= charges[chosen] <= charge_max
    result = []
    for i in chosen[mask]:
        result.append(
            {
                "ion": str(names[i]),
                "moq_u_per_e": float(moqs[i]),
                "mass_u": float(moqs[i] * charges[i]),
                "Z": int(rows[i][4]),
                "offset_percent": float(100 * (moqs[i] / target - 1)),
                **mass_metadata(rows[i], "ame_search"),
            }
        )
    return result
