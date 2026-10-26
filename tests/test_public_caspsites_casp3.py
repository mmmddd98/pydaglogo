from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "examples" / "validate_caspsites_casp3.py"
TEST_OUTPUT_ROOT = ROOT / ".test_runs"


def _load_module():
    spec = importlib.util.spec_from_file_location("validate_caspsites_casp3", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _assert_exact_summary(module, summary: pd.DataFrame, test_type: str) -> None:
    expected = pd.read_csv(module.expected_summary_path_for_test_type(test_type), sep="\t")
    observed = summary[["table", "exact_allclose"]].copy()
    merged = expected.merge(observed, on="table", how="left")

    assert merged["exact_allclose"].notna().all()
    assert (merged["expected_exact_allclose"].astype(bool) == merged["exact_allclose"].astype(bool)).all()
    assert summary["exact_allclose"].all()


@contextmanager
def _managed_output_dir(stem: str):
    TEST_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    output_dir = TEST_OUTPUT_ROOT / f"{stem}_{uuid4().hex}"
    output_dir.mkdir(parents=True, exist_ok=False)
    try:
        yield output_dir
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)


def test_public_validation_api_runs_exact_for_ztest():
    module = _load_module()
    with _managed_output_dir("ztest_api") as output_dir:
        summary, outputs = module.run_validation(output_dir=output_dir, test_type="ztest")

        _assert_exact_summary(module, summary, "ztest")
        assert outputs["summary"].exists()
        assert outputs["diff_logo"].exists()
        assert outputs["motif_logo"].exists()


def test_public_validation_api_runs_exact_for_fisher(exact_fisher_env):
    module = _load_module()
    with _managed_output_dir("fisher_api") as output_dir:
        summary, outputs = module.run_validation(output_dir=output_dir, test_type="fisher")

        _assert_exact_summary(module, summary, "fisher")
        assert outputs["summary"].exists()
        assert outputs["diff_logo"].exists()
        assert outputs["motif_logo"].exists()


def test_public_validation_cli_smoke_ztest():
    with _managed_output_dir("ztest_cli") as output_dir:
        env = os.environ.copy()
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--test-type",
                "ztest",
                "--output-dir",
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )

        assert result.returncode == 0, result.stdout + "\n" + result.stderr
        assert (output_dir / "caspsites_casp3_ztest_parity_summary.tsv").exists()
        assert (output_dir / "caspsites_casp3_ztest_diff_logo.png").exists()
        assert (output_dir / "caspsites_casp3_ztest_motif_logo.png").exists()


def test_public_validation_cli_smoke_fisher(exact_fisher_env):
    with _managed_output_dir("fisher_cli") as output_dir:
        env = os.environ.copy()
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--test-type",
                "fisher",
                "--output-dir",
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )

        assert result.returncode == 0, result.stdout + "\n" + result.stderr
        assert (output_dir / "caspsites_casp3_fisher_parity_summary.tsv").exists()
        assert (output_dir / "caspsites_casp3_fisher_diff_logo.png").exists()
        assert (output_dir / "caspsites_casp3_fisher_motif_logo.png").exists()
