"""The pipeline output must not depend on Python's per-process hash randomisation (PYTHONHASHSEED)."""
import hashlib
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run(out_dir, hash_seed):
    subprocess.run(
        [sys.executable, "-m", "morpho_gruson.cli", "run", "--raw", str(ROOT / "data/raw"),
         "--config", str(ROOT / "config/pipeline.yaml"), "--out", str(out_dir)],
        check=True, capture_output=True, env=dict(os.environ, PYTHONHASHSEED=str(hash_seed)),
    )


def _digests(directory):
    return {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


def test_output_is_identical_for_different_hash_seeds(tmp_path):
    _run(tmp_path / "seed1", 1)
    _run(tmp_path / "seed2", 2)

    first, second = _digests(tmp_path / "seed1"), _digests(tmp_path / "seed2")

    assert first and first == second
