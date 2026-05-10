import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_node_meta_stats_suite_passes():
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not installed in this environment")
    test_files = sorted(str(path.relative_to(ROOT)) for path in (ROOT / "test").glob("*.js"))
    assert test_files
    subprocess.run([node, "--test", *test_files], cwd=ROOT, check=True)


def test_core_module_exports_expected_statistics():
    text = (ROOT / "meta-stats-core.js").read_text(encoding="utf-8")
    required = [
        "export function metaAnalysis",
        "estimators.REML",
        "estimators.DL",
        "estimators.PM",
        "export function radialEgger",
        "export function petersTest",
        "export function petPeese",
        "export function trimFill",
    ]
    missing = [marker for marker in required if marker not in text]
    assert missing == []
