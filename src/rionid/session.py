"""Portable settings, input fingerprints, tables and verified run replay."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .bundles import create_run_directory
from .candidates import matches_selector, parse_ion, split_ions, validate_selectors
from .core import ImportData
from .gui.controller import import_controller
from .io import write_arrays_to_ods
from .version import __version__

SCHEMA_VERSION = 1
MODES = ("Frequency", "Brho", "Kinetic Energy", "Gamma")


@dataclass
class RunConfig:
    datafile: str
    filep: str
    alphap: float
    refion: str
    harmonics: tuple[float, ...]
    mode: str
    value: float
    circumference: float = 108.36
    highlight_ions: str = ""
    reload_data: bool = False
    io_params: dict[str, str] = field(default_factory=dict)
    sim_scalingfactor: Optional[float] = None
    nions: Optional[int] = None
    correct: Optional[tuple[float, ...]] = None
    table_filename: str = "simulation_result.out"
    manual_ions: str = ""
    search_ion: str = ""
    moq_tolerance_percent: float = 1.0
    charge_min: int = 1
    charge_max: Optional[int] = None
    display_filter: str = ""
    frequency_scale_hz_per_100kev: Optional[float] = None
    uncertainty_harmonic: float = 1.0
    show_uncertainties: bool = False
    show_isomers: bool = False
    isomer_filter: str = ""

    def validate(self) -> None:
        if not self.datafile or not self.refion.strip():
            raise ValueError("Select a spectrum and a reference ion")
        if not self.filep and not (self.manual_ions.strip() or self.search_ion.strip()):
            raise ValueError("Select a LISE++ file or provide manual/search candidates")
        parse_ion(self.refion)
        for ion in split_ions(self.manual_ions):
            parse_ion(ion)
        if self.search_ion:
            parse_ion(self.search_ion)
        validate_selectors(self.display_filter)
        validate_selectors(self.isomer_filter)
        validate_selectors(self.highlight_ions)
        if (
            not math.isfinite(self.moq_tolerance_percent)
            or not 0 <= self.moq_tolerance_percent < 100
        ):
            raise ValueError("Mass-to-charge tolerance must be in [0,100) percent")
        if (
            not math.isfinite(self.charge_min)
            or int(self.charge_min) != self.charge_min
            or self.charge_min < 1
        ):
            raise ValueError("Minimum charge must be a positive integer")
        if self.charge_max is not None and (
            not math.isfinite(self.charge_max)
            or int(self.charge_max) != self.charge_max
            or self.charge_max < self.charge_min
        ):
            raise ValueError("Maximum charge must be an integer >= minimum charge")
        if self.frequency_scale_hz_per_100kev is not None and (
            not math.isfinite(self.frequency_scale_hz_per_100kev)
            or self.frequency_scale_hz_per_100kev <= 0
        ):
            raise ValueError("Frequency scale must be finite and positive")
        if not math.isfinite(self.uncertainty_harmonic) or self.uncertainty_harmonic <= 0:
            raise ValueError("Calibration harmonic/multiplier must be finite and positive")
        if self.mode not in MODES:
            raise ValueError(f"Unsupported reference mode: {self.mode}")
        for name in ("alphap", "value", "circumference"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and greater than zero")
        if self.mode == "Gamma" and self.value <= 1:
            raise ValueError("Gamma must be greater than one")
        if not self.harmonics or any(not math.isfinite(h) or h <= 0 for h in self.harmonics):
            raise ValueError("Provide positive finite harmonic orders")
        if len(set(self.harmonics)) != len(self.harmonics):
            raise ValueError("Harmonic orders must not be repeated")
        if self.nions is not None and self.nions <= 0:
            raise ValueError("N ions must be a positive integer")
        if self.sim_scalingfactor is not None and (
            not math.isfinite(self.sim_scalingfactor) or self.sim_scalingfactor <= 0
        ):
            raise ValueError("Scaling factor must be positive and finite")
        if self.correct is not None and (
            len(self.correct) != 3 or any(not math.isfinite(c) for c in self.correct)
        ):
            raise ValueError("Correction requires three finite coefficients: A B C")
        if (
            not self.table_filename
            or Path(self.table_filename).name != self.table_filename
            or self.table_filename
            in (
                ".",
                "..",
                "config.json",
                "candidates.csv",
                "candidates.ods",
                "skipped_candidates.csv",
            )
        ):
            raise ValueError("Text table name must be a filename within the saved-run folder")


def configuration_values(config: RunConfig) -> dict[str, Any]:
    config.validate()
    values = asdict(config)
    for key in ("datafile", "filep"):
        values[key] = str(Path(values[key]).expanduser().resolve()) if values[key] else ""
    values["harmonics"] = list(config.harmonics)
    values["correct"] = list(config.correct) if config.correct is not None else None
    return values


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_inputs(config: RunConfig) -> dict[str, Any]:
    """Identify the actual spectrum/cache and key mapping used by this run."""
    values = configuration_values(config)
    spectrum = Path(values["datafile"])
    cache = spectrum.with_name(spectrum.stem + "_cache.npz")
    cached = not config.reload_data and cache.is_file()
    if cached:
        spectrum = cache
        keys = {"frequency_key": "frequency", "amplitude_key": "amplitude_avg"}
    else:
        keys = {"frequency_key": "arr_0", "amplitude_key": "arr_1", **config.io_params}
    if spectrum.suffix.lower() != ".npz":
        raise ValueError("Saved-run replay currently requires an NPZ spectrum or NPZ cache")
    result = {
        "spectrum": {"path": str(spectrum), "sha256": sha256_file(spectrum), "io_params": keys},
    }
    if values["filep"]:
        result["candidates"] = {"path": values["filep"], "sha256": sha256_file(values["filep"])}
    return result


def validate_inputs(inputs: dict[str, Any]) -> None:
    for name, entry in inputs.items():
        if sha256_file(entry["path"]) != entry["sha256"]:
            raise ValueError(f"Input hash mismatch for {name}: {entry['path']}")


@dataclass
class CompletedRun:
    config: RunConfig
    data: ImportData
    inputs: dict[str, Any]


def run_configuration(
    config: RunConfig,
    *,
    inputs: Optional[dict[str, Any]] = None,
    ionic_inputs: Optional[dict[str, Any]] = None,
) -> CompletedRun:
    """Run without overwriting source caches or exporting into the working directory."""
    config.validate()
    # Retain the settings actually used even if a caller later edits its dataclass.
    config = deepcopy(config)
    inputs = deepcopy(inputs) if inputs is not None else resolve_inputs(config)
    validate_inputs(inputs)
    values = configuration_values(config)
    values.pop("table_filename")
    values.update(
        datafile=inputs["spectrum"]["path"],
        filep=inputs.get("candidates", {}).get("path", ""),
        io_params=inputs["spectrum"]["io_params"],
        reload_data=True,
        harmonics=list(config.harmonics),
        correct=list(config.correct) if config.correct is not None else None,
        output_path=None,
        write_cache=False,
        ionic_inputs=ionic_inputs,
    )
    data = import_controller(**values)
    validate_inputs(inputs)
    return CompletedRun(config, data, inputs)


def frozen_ionic_inputs(data: ImportData) -> dict[str, Any]:
    """Resolved ionic masses preserve the exact inputs used for the overlay."""
    return {
        ion: {
            "moq_u_per_e": float(moq),
            "mass_u": float(data.total_mass[ion]),
            "Z": int(data.protons[ion]),
            "metadata": deepcopy(getattr(data, "candidate_metadata", {}).get(ion, {})),
        }
        for ion, moq in data.moq.items()
    }


def candidate_rows(data: ImportData) -> list[dict[str, Any]]:
    """Export every calculated candidate at each harmonic, including hidden top-N rows."""
    rows = []
    for harmonic in data.simulated_data_dict:
        h = float(harmonic)
        visible = {str(row[2]) for row in data.simulated_data_dict[harmonic]}
        display_filter = getattr(data, "display_filter", "")
        if display_filter:
            visible = {
                name for name in data.nuclei_names if matches_selector(str(name), display_filter)
            }
        for i in np.argsort(data.srrf):
            ion = str(data.nuclei_names[i])
            frev = float(data.srrf[i] * data.ref_frequency)
            metadata = getattr(data, "candidate_metadata", {}).get(ion, {})
            sigma = getattr(data, "frequency_uncertainties", {}).get(harmonic, {}).get(ion)
            rows.append(
                {
                    "ion": ion,
                    "harmonic": h,
                    "frequency_Hz": frev * h,
                    "revolution_frequency_Hz": frev,
                    "relative_frequency": float(data.srrf[i]),
                    "relative_yield": float(data.yield_data[i])
                    if metadata.get("relative_yield_known", True)
                    else None,
                    "source": metadata.get("source", "lise"),
                    "sources": ",".join(metadata.get("sources", ["lise"])),
                    "mass_uncertainty_keV": metadata.get("mass_uncertainty_keV"),
                    "mass_extrapolated": metadata.get("mass_extrapolated"),
                    "frequency_sigma_Hz": sigma,
                    "uncertainty_method": "empirical_mass_scale"
                    if getattr(data, "frequency_scale_hz_per_100kev", None) is not None
                    else "model_mass_derivative",
                    "moq_u_per_e": float(data.moq[ion]),
                    "ionic_mass_u": float(data.total_mass[ion]),
                    "displayed": ion in visible,
                }
            )
    return rows


def verification_values(data: ImportData) -> dict[str, Any]:
    return {
        "reference_frequency_Hz": float(data.ref_frequency),
        "ionic_inputs": frozen_ionic_inputs(data),
        "candidates": candidate_rows(data),
    }


def save_configuration(
    path: str | Path, config: RunConfig, completed: Optional[CompletedRun] = None
) -> None:
    values = configuration_values(config)
    if completed is not None and values != configuration_values(completed.config):
        raise ValueError("These results were computed with different settings")
    inputs = completed.inputs if completed is not None else resolve_inputs(config)
    validate_inputs(inputs)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "software": "RionID",
        "software_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "numpy_version": np.__version__,
        "configuration": values,
        "inputs": inputs,
    }
    if completed is not None:
        payload["ionic_inputs"] = frozen_ionic_inputs(completed.data)
        payload["verification"] = verification_values(completed.data)
        payload["lise_import"] = {
            "summary": getattr(completed.data, "lise_import_summary", {}),
            "skipped_candidates": getattr(completed.data, "skipped_lise_candidates", []),
        }
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_configuration(path: str | Path) -> tuple[RunConfig, dict[str, Any]]:
    source = Path(path).expanduser().resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    if payload.get("schema_version") != SCHEMA_VERSION or payload.get("software") != "RionID":
        raise ValueError("Unsupported RionID configuration schema or software")
    values = dict(payload["configuration"])
    for key in ("datafile", "filep"):
        if values[key]:
            p = Path(values[key]).expanduser()
            values[key] = str(p if p.is_absolute() else source.parent / p)
    values["harmonics"] = tuple(values["harmonics"])
    if values.get("correct") is not None:
        values["correct"] = tuple(values["correct"])
    config = RunConfig(**values)
    config.validate()
    for entry in payload["inputs"].values():
        p = Path(entry["path"]).expanduser()
        entry["path"] = str(p if p.is_absolute() else source.parent / p)
    return config, payload


def _compare(expected: Any, actual: Any, label: str = "results") -> list[str]:
    if isinstance(expected, dict) and isinstance(actual, dict):
        if not expected.keys() <= actual.keys():
            return [f"{label}: different fields"]
        return [msg for k in expected for msg in _compare(expected[k], actual[k], f"{label}.{k}")]
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            return [f"{label}: different row count"]
        return [
            msg
            for i, (a, b) in enumerate(zip(expected, actual))
            for msg in _compare(a, b, f"{label}[{i}]")
        ]
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        if isinstance(actual, (int, float)) and not isinstance(actual, bool):
            atol = 1e-5 if label.endswith("_Hz") else 1e-12
            if math.isclose(
                expected, actual, rel_tol=0 if label.endswith("_Hz") else 1e-12, abs_tol=atol
            ):
                return []
    elif expected == actual:
        return []
    return [f"{label}: expected {expected!r}, obtained {actual!r}"]


@dataclass
class ReplayReport:
    completed: CompletedRun
    warnings: list[str]
    mismatches: list[str]

    @property
    def verified(self) -> bool:
        return not self.warnings and not self.mismatches


def replay_configuration(
    path: str | Path, *, input_overrides: Optional[dict[str, str]] = None
) -> ReplayReport:
    config, payload = load_configuration(path)
    for name, replacement in (input_overrides or {}).items():
        if name not in payload["inputs"]:
            raise ValueError(f"Unknown input: {name}")
        payload["inputs"][name]["path"] = str(Path(replacement).expanduser().resolve())
        if name == "candidates":
            config.filep = payload["inputs"][name]["path"]
        elif name == "spectrum":
            config.datafile = payload["inputs"][name]["path"]
            config.io_params = dict(payload["inputs"][name]["io_params"])
            config.reload_data = True
    completed = run_configuration(
        config, inputs=payload["inputs"], ionic_inputs=payload.get("ionic_inputs")
    )
    if "lise_import" in payload:
        completed.data.lise_import_summary = payload["lise_import"]["summary"]
        completed.data.skipped_lise_candidates = payload["lise_import"]["skipped_candidates"]
    warnings = []
    if payload.get("software_version") != __version__:
        warnings.append("Software version differs from the saved run")
    if payload.get("numpy_version") != np.__version__:
        warnings.append("NumPy version differs from the saved run")
    mismatches = []
    if "verification" in payload:
        mismatches = _compare(payload["verification"], verification_values(completed.data))
    else:
        warnings.append("Settings-only configuration: no numerical checkpoints to verify")
    return ReplayReport(completed, warnings, mismatches)


def export_tables(completed: CompletedRun, out: Path) -> list[Path]:
    data = completed.data
    rows = candidate_rows(data)
    columns = list(rows[0])
    csv_path = out / "candidates.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    ods_path = out / "candidates.ods"
    write_arrays_to_ods(
        str(ods_path),
        "Candidates",
        columns,
        *[[row[c] if row[c] is not None else "" for row in rows] for c in columns],
    )
    text_path = out / completed.config.table_filename
    with text_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    written = [csv_path, ods_path, text_path]
    skipped = getattr(data, "skipped_lise_candidates", [])
    if skipped:
        skipped_path = out / "skipped_candidates.csv"
        with skipped_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["record", "ion", "kind", "reason"])
            writer.writeheader()
            writer.writerows(skipped)
        written.append(skipped_path)
    return written


def save_bundle(
    parent: str | Path,
    label: str,
    config: RunConfig,
    completed: Optional[CompletedRun] = None,
) -> tuple[Path, list[Path]]:
    # Validate before allocating a folder, including stale source files.
    values = configuration_values(config)
    if completed is not None and values != configuration_values(completed.config):
        raise ValueError("These results were computed with different settings")
    validate_inputs(completed.inputs if completed is not None else resolve_inputs(config))
    out = create_run_directory(parent, "RionID", label)
    config_path = out / "config.json"
    save_configuration(config_path, config, completed)
    written = [config_path]
    if completed is not None:
        written.extend(export_tables(completed, out))
    return out, written
