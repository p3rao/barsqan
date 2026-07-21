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
"""
from __future__ import annotations

import csv
import glob
import os
from collections import defaultdict
from typing import Dict, List


def load_all_counts(csv_dir: str, min_count: int = 2) -> Dict[str, Dict[str, int]]:
    """Returns {barcode: {sample: count}} across all *.csv files in csv_dir."""
    barcode_sample_counts: Dict[str, Dict[str, int]] = defaultdict(dict)

    for path in sorted(glob.glob(os.path.join(csv_dir, "*.csv"))):
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample = row["sample"]
                bc = row["corrected_barcode"]
                count = int(row["total_count"])

                if "N" in bc.upper():
                    continue
                if count < min_count:
                    continue

                barcode_sample_counts[bc][sample] = barcode_sample_counts[bc].get(sample, 0) + count

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
