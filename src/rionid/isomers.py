"""Nubase2020 isomers: excited states listed with an excitation energy."""

import os
import re
import urllib.request as ur
from dataclasses import dataclass

from .masses import AMEData

ISOMER_LABELS = "mnpqrx"


@dataclass(frozen=True)
class Isomer:
    label: str
    excitation_keV: float
    excitation_unc_keV: float
    extrapolated: bool


class NubaseData:
    """Loads and caches the Nubase2020 isomer table.

    On construction, reads `~/.ame/nubase.data`, downloading it first if
    absent. Only isomers are kept (Nubase states m, n and, for nuclides with
    more than two, p, q, r, x), indexed by (element name, mass number).
    Nubase level (p, q) and resonance rows carry no half-life and are skipped.
    """

    NUBASE_DATA_LINK = "https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt"

    def __init__(self, path=None):
        folder = os.path.expanduser("~") + AMEData.FOLDER_NAME
        os.makedirs(folder, exist_ok=True)
        self.nubase_data_filename = path or f"{folder}nubase.data"
        if not os.path.exists(self.nubase_data_filename):
            self._download()
        self._isomers = {}
        self._parse_nubase()

    def isomers(self, name, aa):
        """Returns the isomers of (element name, mass number `aa`), lowest label first."""
        return self._isomers.get((name, aa), [])

    def _download(self):
        req = ur.Request(self.NUBASE_DATA_LINK, headers={"User-Agent": "Magic Browser"})
        with open(self.nubase_data_filename, "b+w") as f:
            f.write(ur.urlopen(req).read())

    def _parse_nubase(self):
        """Parses the fixed-width Nubase2020 table; header lines are skipped."""
        with open(self.nubase_data_filename) as f:
            for line in f:
                if not (line[:3].isdigit() and line[4:8].isdigit()):
                    continue
                state, label = int(line[7]), line[16]
                if not 1 <= state <= 6 or label not in ISOMER_LABELS:
                    continue
                if state > 2 and not line[69:78].strip():
                    continue
                energy, sigma = line[42:54].strip(), line[54:65].strip()
                try:
                    excitation = float(energy.replace("#", ""))
                    uncertainty = float(sigma.replace("#", "") or 0)
                except ValueError:
                    continue
                if excitation <= 0:
                    continue
                element = re.sub(r"^\d+", "", line[11:16].strip())
                self._isomers.setdefault((element, int(line[:3])), []).append(
                    Isomer(label, excitation, uncertainty, "#" in energy + sigma)
                )


_nubase_cache = None


def get_nubase_data():
    """Returns a process-lifetime-cached `NubaseData` instance."""
    global _nubase_cache
    if _nubase_cache is None:
        _nubase_cache = NubaseData()
    return _nubase_cache
