"""Read LISE++ candidates and report unavailable masses once per input file."""

import logging

import numpy as np

from rionid.candidates import parse_ion
from rionid.masses import get_ame_data

log = logging.getLogger(__name__)
_reported = set()


class MissingMassError(ValueError):
    """A well-formed candidate has no mass in the selected AME table."""


class LISEreader:
    def __init__(self, filename, *, validate_masses=True):
        self.ame = get_ame_data() if validate_masses else None
        self.ame_data = self.ame.ame_table if self.ame is not None else []
        self.skipped_entries = []
        self.summary = {}
        self.filename = str(filename)
        self._read(filename)
        self._indices = {}
        for i, row in enumerate(self.data):
            if len(row) < 2:
                continue
            try:
                name = parse_ion(row[0] + "+" + row[1]).name
            except ValueError:
                continue
            self._indices.setdefault(name, i)

    def _read(self, filename):
        with open(filename, encoding="latin1") as f:
            lines = f.readlines()

        sections = [i + 1 for i, line in enumerate(lines) if "[Calculations]" in line]
        if not sections or sections[-1] >= len(lines):
            raise ValueError("LISE++ file has no nonempty [Calculations] section")
        calculation_lines = [
            line
            for line in lines[sections[-1] :]
            if line.strip() and not line.lstrip().startswith(("#", ";"))
        ]
        first_record = None
        for line in calculation_lines:
            fields = line.replace("+", "").split()
            if len(fields) < 3:
                continue
            try:
                parse_ion(fields[0] + "+" + fields[1])
            except ValueError:
                continue
            first_record = line
            break
        if first_record is None:
            raise ValueError("LISE++ file has no parsable calculation records")
        self.centre_index = len(first_record.split()) - 1

        self.data = []
        for line in calculation_lines:
            fields = line.replace("+", "").split()
            if len(fields) <= self.centre_index:
                # Keep damaged rows so the diagnostic list accounts for them.
                self.data.append(fields)
                continue
            values = line.replace("=", "").split()[self.centre_index].split(",")
            self.data.append(fields[: self.centre_index] + values)

    def get_index(self, name):
        canonical = parse_ion(name).name
        try:
            return self._indices[canonical]
        except KeyError:
            raise ValueError(f"{canonical} is not present in the LISE++ file") from None

    def get_info(self, name):
        ion = parse_ion(name)
        index = self.get_index(name)
        record = self.data[index]
        if len(record) <= self.centre_index:
            raise ValueError(f"{ion.name} has an incomplete transport or yield record")
        charges = list(map(int, record[1 : self.centre_index]))
        yield_value = float(record[self.centre_index])
        # Intermediate transport columns can retain the original beam setup.
        # Only the final charge defines the ion whose mass is simulated.
        if not charges or not 1 <= charges[-1] <= ion.Z:
            raise ValueError(f"{ion.name} has an invalid final charge state in the LISE++ file")
        if not np.isfinite(yield_value) or yield_value < 0:
            raise ValueError(f"{ion.name} has an invalid LISE++ yield")
        if self.ame is not None:
            row = self.ame.lookup(ion.element, ion.A)
            if row is None:
                raise MissingMassError(f"No AME2020 mass is available for {ion.A}{ion.element}")
            nuclide = [row[6], row[5], row[4], row[3]]
        else:
            # Frozen replay already contains the ionic masses; no live AME lookup is needed.
            nuclide = [ion.element, ion.A, ion.Z, ion.A - ion.Z]
        return nuclide + [charges, yield_value]

    def get_lise_all(self):
        data = np.array(self.data)
        return np.transpose([data[:, 0], data[:, 5], data[:, 6]])

    def get_info_all(self):
        """Keep the legacy first-match policy for duplicate transport records."""
        result = []
        self.skipped_entries = []
        missing = set()
        for i, line in enumerate(self.data):
            name = line[0] + "+" + line[1] if len(line) >= 2 else " ".join(line)
            try:
                result.append(self.get_info(name))
            except ValueError as error:
                kind = "missing_mass" if isinstance(error, MissingMassError) else "invalid_record"
                if kind == "missing_mass":
                    ion = parse_ion(name)
                    missing.add(f"{ion.A}{ion.element}")
                self.skipped_entries.append(
                    {
                        "record": i + 1,
                        "ion": name,
                        "kind": kind,
                        "reason": str(error),
                    }
                )
        self.summary = {
            "input_records": len(self.data),
            "usable_records": len(result),
            "skipped_records": len(self.skipped_entries),
            "missing_mass_nuclides": sorted(missing),
            "invalid_records": sum(r["kind"] == "invalid_record" for r in self.skipped_entries),
        }
        # A run re-reads the same file (preview, replay); warn once per file and result.
        key = (self.filename, len(result), len(self.skipped_entries))
        if self.skipped_entries:
            report = log.debug if key in _reported else log.warning
            _reported.add(key)
            report(
                "LISE++: %d usable records; omitted %d records (%d nuclides absent from AME2020, "
                "%d invalid records). Save results for the detailed skipped-candidate list.",
                len(result),
                len(self.skipped_entries),
                len(missing),
                self.summary["invalid_records"],
            )
        return result

    def get_info_specific(self, param_index_list):
        return_list = [
            [LISEreader.float_check(line[i]) for i in param_index_list] for line in self.data
        ]
        return return_list

    def get_yield_sorted(self):
        name_yield = [
            [str(line[0]) + "+" + str(line[1]), LISEreader.float_check(line[7])]
            for line in self.data
        ]
        sorted_by_yield = sorted(name_yield, key=lambda x: x[1])[::-1]
        return sorted_by_yield

    @staticmethod
    def float_check(value):  # returns float if not string
        if value.replace(".", "").replace("e-", "").replace("e+", "").isdigit():
            return float(value)
        else:
            return value
