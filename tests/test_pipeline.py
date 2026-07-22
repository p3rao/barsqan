"""
Minimal smoke test for the full pipeline, using the synthetic FASTQ data in
example_data/. Run with:  python3 -m pytest tests/  (or just `python3 tests/test_pipeline.py`)
"""
import csv
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE_DIR = os.path.join(ROOT, "example_data")


def run(cmd):
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed: {cmd}\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}"
    return result


def test_full_pipeline(tmp_path=None):
    outdir = os.path.join(EXAMPLE_DIR, "_test_results")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    run([
        sys.executable, "-m", "barsqan.cli", "run",
        "--samples", os.path.join(EXAMPLE_DIR, "samples.txt"),
        "--fastq-dir", EXAMPLE_DIR,
        "--outdir", outdir,
        "--index-map", os.path.join(EXAMPLE_DIR, "index_map.csv"),
        "--do-map",
        "--config", os.path.join(ROOT, "config.example.yaml"),
    ])

    # extraction: Sample_A has 15 reads engineered to fail (5 each stage) out of 70
    stats_path = os.path.join(outdir, "1_extracted", "Sample_A_S1_L001.stats.tsv")
    stats = dict(line.rstrip("\n").split("\t") for line in open(stats_path) if line.startswith("#"))
    assert stats["# Total_reads"] == "70"
    assert stats["# Passed"] == "55"
    assert stats["# Quality_fail"] == "5"
    assert stats["# Regex_fail"] == "5"
    assert stats["# Overlap_fail"] == "5"

    # clustering: Sample_A should collapse to exactly 3 distinct barcodes
    with open(os.path.join(outdir, "3_clustered", "Sample_A_S1_L001.corrected_barcodes.csv")) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    counts = {r["corrected_barcode"]: int(r["total_count"]) for r in rows}
    assert sorted(counts.values(), reverse=True) == [30, 15, 10]

    # mapping: exactly 2 barcodes shared between Sample_A and Sample_B
    with open(os.path.join(outdir, "4_mapped", "barcode_overlap_summary.csv")) as f:
        map_rows = list(csv.DictReader(f))
    shared = [r for r in map_rows if int(r["n_samples"]) == 2]
    assert len(shared) == 2

    shutil.rmtree(outdir)
    print("OK: full pipeline smoke test passed")


if __name__ == "__main__":
    test_full_pipeline()
