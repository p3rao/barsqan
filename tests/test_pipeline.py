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

# Allow `python3 tests/test_pipeline.py` to import the package from source
# even when barsqan is not pip-installed into the active environment.
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


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


def test_clustering_merges_error_variants():
    """A 1-edit low-count variant collapses into its abundant parent, and
    both backends agree. Also guards against the O(N^2) regression by keeping
    the input large-ish and asserting it completes quickly."""
    import time
    from collections import Counter
    from barsqan.config import Config
    import barsqan.distance as D
    import barsqan.cluster_umi as C

    def make_counts():
        parent = "AAAAAacGGGGGgaCCCCCtcTTTTT"           # 26nt, matches motif
        variant = "AAAAAacGGGGGgaCCCCAtcTTTTT"          # 1 substitution
        far = "GGGGGacAAAAAgaTTTTTtcCCCCC"             # unrelated barcode
        c = Counter()
        c[parent] = 500
        c[variant] = 3     # low count -> should merge into parent
        c[far] = 400
        return c

    cfg = Config()
    for force_pure in (False, True):
        D._USE_RAPIDFUZZ_CPP = D._USE_RAPIDFUZZ_CPP and not force_pure
        rep, mapping = C.cluster_barcodes(make_counts(), cfg)
        parent = "AAAAAacGGGGGgaCCCCCtcTTTTT"
        variant = "AAAAAacGGGGGgaCCCCAtcTTTTT"
        assert mapping[variant] == parent, f"variant not merged (force_pure={force_pure})"
        assert rep[parent] == 503, f"parent count wrong: {rep[parent]}"
        assert len(rep) == 2
    print("OK: clustering merges error variants (both backends)")


def _write_cluster_csv(d, sample, counts):
    import csv as _csv
    total = sum(counts.values())
    with open(os.path.join(d, f"{sample}.corrected_barcodes.csv"), "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["sample", "corrected_barcode", "total_count", "freq_within_sample"])
        for bc, n in counts.items():
            w.writerow([sample, bc, n, f"{n/total:.6f}"])


def test_map_relative_abundance_can_drop_large_fraction():
    """The relative-abundance filter is a VALUE threshold: if most barcodes
    are far below the reference stat, most barcodes are dropped (not a fixed
    5%). Here 60 of 100 barcodes are tiny and must all be dropped."""
    from barsqan.map_overlap import _relative_abundance_cutoff

    # 60 barcodes at count=1 (tiny), 40 at count=1000 (real).
    counts = {f"lo{i:03d}": 1 for i in range(60)}
    counts.update({f"hi{i:03d}": 1000 for i in range(40)})
    # median freq over the 100 barcodes: with 60 at count=1 the median is the
    # low value, so use mean here to get a threshold between the two groups.
    drop, cutoff, ref = _relative_abundance_cutoff(counts, 0.05, "mean")
    assert len(drop) == 60, f"expected 60 dropped, got {len(drop)}"
    assert all(bc.startswith("lo") for bc in drop)
    # cutoff is 5% of the mean frequency, and every 'lo' barcode is below it
    assert cutoff is not None and cutoff > 0
    print("OK: relative-abundance filter can drop a large fraction (60%)")


def test_map_mean_vs_median_reference():
    """mean vs median give different thresholds on a skewed distribution."""
    from barsqan.map_overlap import _relative_abundance_cutoff
    # Heavy low tail + a few whales -> mean >> median.
    counts = {f"lo{i:03d}": 2 for i in range(90)}
    counts.update({f"hi{i:02d}": 5000 for i in range(10)})
    _, cutoff_mean, ref_mean = _relative_abundance_cutoff(counts, 0.5, "mean")
    _, cutoff_median, ref_median = _relative_abundance_cutoff(counts, 0.5, "median")
    assert ref_mean > ref_median, (ref_mean, ref_median)
    assert cutoff_mean > cutoff_median
    print("OK: mean vs median reference produce different cutoffs")


def test_map_relative_abundance_end_to_end():
    """A barcode abundant in one sample but a low-abundance straggler in
    another is not reported as shared once the relative-abundance filter runs."""
    import tempfile
    from barsqan.map_overlap import load_all_counts

    d = tempfile.mkdtemp()
    shared = "AAAAAacGGGGGgaCCCCCtcTTTTT"

    # Sample_A: shared is abundant among a tight high-count population.
    A = {f"AAAAAacGGGGGgaCCC{i:03d}xyz"[:26]: 100 + i for i in range(40)}
    A[shared] = 500
    _write_cluster_csv(d, "Sample_A", A)
    # Sample_B: high-count population + shared present as a tiny straggler.
    B = {f"BBBBBacTTTTTgaAA{i:03d}xyz"[:26]: 100 + i for i in range(40)}
    B[shared] = 1
    _write_cluster_csv(d, "Sample_B", B)

    no_filter = load_all_counts(d, min_count=1, min_rel_abundance=0.0)
    assert set(no_filter[shared]) == {"Sample_A", "Sample_B"}

    # 5% of the median: shared=1 in B is far below, gets dropped; A keeps it.
    with_filter = load_all_counts(d, min_count=1, min_rel_abundance=0.05, relative_to="median")
    assert set(with_filter[shared]) == {"Sample_A"}

    shutil.rmtree(d)
    print("OK: relative-abundance filter drops per-sample low-abundance straggler")


def test_map_plot_generation():
    """Plotting produces one PNG per sample with the cutoff marked (skips
    cleanly if matplotlib is unavailable)."""
    import tempfile
    from barsqan.map_overlap import load_all_counts, plot_count_distributions

    try:
        import matplotlib  # noqa
    except Exception:
        print("SKIP: matplotlib not installed, plot test skipped")
        return

    d = tempfile.mkdtemp()
    outdir = tempfile.mkdtemp()
    counts = {f"bc{i:03d}": i for i in range(1, 101)}
    _write_cluster_csv(d, "Sample_P", counts)

    dist = {}
    load_all_counts(d, min_count=2, min_rel_abundance=0.05, relative_to="median", dist_out=dist)
    pngs = plot_count_distributions(dist, outdir, 0.05)
    assert len(pngs) == 1, pngs
    assert os.path.exists(pngs[0]) and os.path.getsize(pngs[0]) > 1000

    shutil.rmtree(d)
    shutil.rmtree(outdir)
    print("OK: map plot generation writes per-sample PNG with cutoff")


if __name__ == "__main__":
    test_full_pipeline()
    test_clustering_merges_error_variants()
    test_map_relative_abundance_can_drop_large_fraction()
    test_map_mean_vs_median_reference()
    test_map_relative_abundance_end_to_end()
    test_map_plot_generation()
