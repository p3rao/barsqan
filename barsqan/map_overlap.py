"""
Step 4 (optional): map barcode occurrences across multiple sample-name
count tables produced by cluster_umi.cluster_sample.

Builds a barcode x sample matrix of counts and within-sample frequencies,
plus a per-barcode summary of how many samples it appears in. This is the
generic replacement for the original Plate/Row/Col-specific overlap script:
it works for arbitrary sample names, and if your samples DO follow a
Plate/Row/Col-style naming convention you can still recover that breakdown
by grouping on a substring of the sample name after the fact (the full
per-sample matrix retains everything needed to do that downstream).

Barcodes containing 'N' are excluded (ambiguous base calls), and a barcode
is only counted for a given sample if its count in that sample is >= 2
(drops singleton reads which are more likely to be errors than real
lineages) - matching the conventions used in the earlier lab scripts.

Optionally, before pooling a sample's barcodes into the cross-sample matrix,
the bottom X%% of that sample's *own* barcode count distribution can be
dropped (`bottom_percentile`). This removes the long low-count tail (likely
sequencing error / contamination) on a per-sample basis, so a barcode that
is abundant in one sample but only a low-count straggler in another is not
spuriously reported as "shared".
"""
from __future__ import annotations

import csv
import glob
import os
from collections import defaultdict
from typing import Dict, List, Tuple


def _percentile_threshold(sorted_counts: List[int], pct: float) -> float:
    """Value of the pct-th percentile of an ascending-sorted list of counts.

    Uses linear interpolation between closest ranks (same convention as
    numpy.percentile with the default 'linear' method), but implemented in
    pure Python so map has no numpy dependency. Barcodes with count <= this
    threshold are considered "bottom pct%" and dropped.
    """
    if not sorted_counts:
        return float("-inf")
    if len(sorted_counts) == 1:
        # a single barcode: never drop it via the percentile filter
        return float("-inf")
    rank = (pct / 100.0) * (len(sorted_counts) - 1)
    lo = int(rank)
    hi = min(lo + 1, len(sorted_counts) - 1)
    frac = rank - lo
    return sorted_counts[lo] + (sorted_counts[hi] - sorted_counts[lo]) * frac


def _read_sample_counts(path: str, min_count: int) -> Tuple[str, Dict[str, int]]:
    """Read one cluster CSV -> (sample_name, {barcode: count}).

    Applies the N-base and absolute min_count filters, but NOT the
    per-sample percentile filter (that is applied afterwards, once the
    whole sample distribution is known).
    """
    sample_name = None
    counts: Dict[str, int] = defaultdict(int)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["sample"]
            sample_name = sample_name or sample
            bc = row["corrected_barcode"]
            count = int(row["total_count"])
            if "N" in bc.upper():
                continue
            if count < min_count:
                continue
            counts[bc] += count
    return sample_name, dict(counts)


def load_all_counts(
    csv_dir: str,
    min_count: int = 2,
    bottom_percentile: float = 0.0,
) -> Dict[str, Dict[str, int]]:
    """Returns {barcode: {sample: count}} across all *.csv files in csv_dir.

    min_count:         absolute minimum per-sample count for a barcode
                       (applied first, per barcode).
    bottom_percentile: if > 0, for each sample independently, drop barcodes
                       whose count is at or below the bottom_percentile-th
                       percentile of that sample's barcode count
                       distribution (applied after min_count, before pooling
                       barcodes across samples). e.g. 5.0 drops the bottom 5%.
    """
    barcode_sample_counts: Dict[str, Dict[str, int]] = defaultdict(dict)

    for path in sorted(glob.glob(os.path.join(csv_dir, "*.csv"))):
        sample_name, counts = _read_sample_counts(path, min_count)
        if not counts:
            continue

        threshold = float("-inf")
        if bottom_percentile and bottom_percentile > 0:
            sorted_counts = sorted(counts.values())
            threshold = _percentile_threshold(sorted_counts, bottom_percentile)
            kept = {bc: c for bc, c in counts.items() if c > threshold}
            dropped = len(counts) - len(kept)
            import sys
            print(
                f"[map] {sample_name}: dropped bottom {bottom_percentile:g}% "
                f"({dropped}/{len(counts)} barcodes with count <= {threshold:.2f})",
                file=sys.stderr,
            )
            counts = kept

        for bc, count in counts.items():
            barcode_sample_counts[bc][sample_name] = (
                barcode_sample_counts[bc].get(sample_name, 0) + count
            )

    return barcode_sample_counts


def write_overlap_outputs(
    barcode_sample_counts: Dict[str, Dict[str, int]],
    outdir: str,
) -> None:
    os.makedirs(outdir, exist_ok=True)

    all_samples: List[str] = sorted({s for counts in barcode_sample_counts.values() for s in counts})
    sample_totals: Dict[str, int] = defaultdict(int)
    for counts in barcode_sample_counts.values():
        for s, c in counts.items():
            sample_totals[s] += c

    counts_path = os.path.join(outdir, "barcode_x_sample_counts.csv")
    freq_path = os.path.join(outdir, "barcode_x_sample_freq.csv")
    summary_path = os.path.join(outdir, "barcode_overlap_summary.csv")

    with open(counts_path, "w", newline="") as f_counts, open(freq_path, "w", newline="") as f_freq:
        w_counts = csv.writer(f_counts)
        w_freq = csv.writer(f_freq)
        w_counts.writerow(["barcode"] + all_samples)
        w_freq.writerow(["barcode"] + all_samples)

        # sort barcodes by number of samples present (descending), then total count
        def sort_key(item):
            bc, counts = item
            return (-len(counts), -sum(counts.values()))

        for bc, counts in sorted(barcode_sample_counts.items(), key=sort_key):
            count_row = [counts.get(s, "") for s in all_samples]
            freq_row = [
                f"{counts[s] / sample_totals[s]:.6f}" if s in counts and sample_totals[s] else ""
                for s in all_samples
            ]
            w_counts.writerow([bc] + count_row)
            w_freq.writerow([bc] + freq_row)

    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["barcode", "n_samples", "samples", "total_count"])
        for bc, counts in sorted(
            barcode_sample_counts.items(),
            key=lambda item: (-len(item[1]), -sum(item[1].values())),
        ):
            w.writerow([bc, len(counts), ";".join(sorted(counts)), sum(counts.values())])
