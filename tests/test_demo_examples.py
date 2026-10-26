from __future__ import annotations

import importlib.util
import shutil
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[1]
ZTEST_SCRIPT = ROOT / "examples" / "demo_flanks_oneshot_ztest.py"
FISHER_SCRIPT = ROOT / "examples" / "demo_flanks_oneshot_fisher.py"
TEST_OUTPUT_ROOT = ROOT / ".test_runs"


def _load_module(name: str, script_path: Path):
    spec = importlib.util.spec_from_file_location(name, script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@contextmanager
def _managed_output_dir(stem: str):
    TEST_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    output_dir = TEST_OUTPUT_ROOT / f"{stem}_{uuid4().hex}"
    output_dir.mkdir(parents=True, exist_ok=False)
    try:
        yield output_dir
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)


def test_flank_demo_ztest_smoke():
    module = _load_module("demo_flanks_oneshot_ztest", ZTEST_SCRIPT)
    with _managed_output_dir("flank_ztest_demo") as output_dir:
        module.OUTPUT_DIR = output_dir

        exit_code = module.main()

        assert exit_code == 0
        assert (module.OUTPUT_DIR / "demo_flanks_oneshot_ztest_diff_logo.png").exists()
        assert (module.OUTPUT_DIR / "demo_flanks_oneshot_ztest_motif_logo.png").exists()


def test_flank_demo_fisher_smoke():
    module = _load_module("demo_flanks_oneshot_fisher", FISHER_SCRIPT)
    with _managed_output_dir("flank_fisher_demo") as output_dir:
        module.OUTPUT_DIR = output_dir

        exit_code = module.main()

        assert exit_code == 0
        assert (module.OUTPUT_DIR / "demo_flanks_oneshot_fisher_diff_logo.png").exists()
        assert (module.OUTPUT_DIR / "demo_flanks_oneshot_fisher_motif_logo.png").exists()
