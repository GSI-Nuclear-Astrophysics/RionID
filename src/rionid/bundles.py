"""Create dated run folders without replacing previous runs."""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path


def create_run_directory(parent: str | Path, software: str, label: str = "") -> Path:
    """Return a newly created DATE_SOFTWARE_CONFIG[_label][_NN] folder."""
    suffix = re.sub(r"[^\w.-]+", "_", label.strip(), flags=re.UNICODE).strip("_.")
    name = f"{date.today().isoformat()}_{software}_CONFIG" + (f"_{suffix}" if suffix else "")
    root = Path(parent).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    index = 1
    while True:
        candidate = root / (name if index == 1 else f"{name}_{index:02d}")
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            index += 1
