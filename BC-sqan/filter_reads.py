"""
Step 2: filter extracted barcodes by internal index + a tolerant re-check
of the barcode motif.

Even though extraction already enforces the barcode motif strictly (0
mismatches allowed in the fixed bases) and a rough index window is implied
by read layout, this step provides an independent, tolerant sanity check
against a *known expected index per sample* (e.g. from a plate/well index
map) and re-verifies the barcode motif with a configurable mismatch
tolerance. This catches cases like index cross-talk / bleed between
samples that the extraction step alone would not detect.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Optional

from .config import Config
from .motif import barcode_matches_motif


def hamming(s1: str, s2: str) -> Optional[int]:
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def index_matches(observed: str, expected: str, chunk_len: int, max_mismatch: int) -> bool:
    """Check observed vs expected index in independent chunks of chunk_len,
    each tolerating up to max_mismatch. Chunking prevents one systematic
    mismatch in one half of a combined index from being masked by a perfect
    match in the other half.
    """
    if len(observed) != len(expected):
        return False
    for start in range(0, len(observed), chunk_len):
        o = observed[start:start + chunk_len]
        e = expected[start:start + chunk_len]
        d = hamming(o, e)
        if d is None or d > max_mismatch:
            return False
    return True


def load_expected_index_map(mapping_csv: str) -> Dict[str, str]:
    """CSV with columns: sample,expected_index"""
    expected = {}
    with open(mapping_csv) as f:
        reader = csv.DictReader(f)
        if "sample" not in reader.fieldnames or "expected_index" not in reader.fieldnames:
            raise ValueError(
                f"{mapping_csv} must have columns 'sample' and 'expected_index' "
                f"(found: {reader.fieldnames})"
            )
        for row in reader:
            expected[row["sample"]] = row["expected_index"].strip().upper()
    return expected


@dataclass
class FilterCounts:
    kept: int = 0
    index_mismatch: int = 0
    motif_absent: int = 0


def filter_sample(
    parsed_tsv: str,
    sample_name: str,
    expected_index: str,
    cfg: Config,
    kept_out_dir: str,
    filtered_out_dir: str,
) -> FilterCounts:
    os.makedirs(kept_out_dir, exist_ok=True)
    os.makedirs(filtered_out_dir, exist_ok=True)

    kept_path = os.path.join(kept_out_dir, f"{sample_name}.kept.tsv")
    filt_path = os.path.join(filtered_out_dir, f"{sample_name}.filtered.tsv")
    counts = FilterCounts()

    with open(parsed_tsv) as fin, open(kept_path, "w") as fout_kept, open(filt_path, "w") as fout_filt:
        header = fin.readline().rstrip("\n")
        cols = header.split("\t")
        try:
            idx_col = cols.index("Combined_Index")
            bc_col = cols.index("Barcode")
        except ValueError:
            raise ValueError(f"Required columns not found in {parsed_tsv} (need Combined_Index, Barcode)")

        new_header = header + "\tfilter_reason"
        fout_kept.write(new_header + "\n")
        fout_filt.write(new_header + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            combined_idx = parts[idx_col].upper()
            barcode = parts[bc_col]

            if not index_matches(combined_idx, expected_index, cfg.index_chunk_len, cfg.index_max_mismatch):
                fout_filt.write(line + "\tindex_mismatch\n")
                counts.index_mismatch += 1
                continue

            if not barcode_matches_motif(barcode, cfg.barcode_motif, cfg.motif_max_mismatch):
                fout_filt.write(line + "\tmotif_absent\n")
                counts.motif_absent += 1
                continue

            fout_kept.write(line + "\t\n")
            counts.kept += 1

    return counts


def write_summaries(
    per_sample_counts: Dict[str, FilterCounts],
    summary_csv: str,
) -> None:
    with open(summary_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "kept", "index_mismatch", "motif_absent", "total", "kept_fraction"])
        for sample, c in per_sample_counts.items():
            total = c.kept + c.index_mismatch + c.motif_absent
            frac = c.kept / total if total else 0.0
            w.writerow([sample, c.kept, c.index_mismatch, c.motif_absent, total, f"{frac:.4f}"])
