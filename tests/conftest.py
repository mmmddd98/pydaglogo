from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pydaglogo.statistics import _candidate_rscript_paths


def _resolved_rscript_path() -> Path | None:
    override = os.environ.get("PYDAGLOGO_RSCRIPT")
    if override:
        path = Path(override)
        if path.exists():
            return path

    candidates = _candidate_rscript_paths()
    if candidates:
        return candidates[0]

    return None


@pytest.fixture
def exact_fisher_env(monkeypatch) -> Path:
    rscript = _resolved_rscript_path()
    if rscript is None:
        pytest.skip("Exact Fisher validation requires a working Rscript executable.")

    monkeypatch.setenv("PYDAGLOGO_FISHER_BACKEND", "r")
    monkeypatch.setenv("PYDAGLOGO_RSCRIPT", str(rscript))
    return rscript
