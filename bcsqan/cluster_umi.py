"""
Step 3: collapse UMIs and cluster barcode sequencing-error variants into a
final per-sample barcode count table.

Two-stage collapsing:
  1. UMI collapsing: within a given barcode, group UMIs that are within
     `umi_group_max_dist` Hamming distance of each other (keeps distinct
     true molecules while absorbing sequencing errors in the UMI).
  2. Barcode clustering: across all (barcode, umi_group) counts, greedily
     assign low-abundance barcode variants that are within
     `cluster_max_dist` Levenshtein distance of a much higher-abundance
     "parent" barcode (count <= cluster_merge_ratio * parent_count) to that
     parent, unless the variant's own count already exceeds
     cluster_protect_count (in which case it is kept as its own barcode).
"""
from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from typing import Dict, Tuple

from rapidfuzz.distance import Levenshtein  # pip install rapidfuzz

from .config import Config


def hamming(s1: str, s2: str):
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def load_barcode_umi_counts(tsv_path: str, cfg: Config) -> Counter:
    """Load barcode + UMI-group counts from a kept/parsed TSV.

    Returns Counter keyed by "{barcode}_{umi_group_repr}" -> read count.
    """
    bc_umi_groups: Dict[str, Counter] = defaultdict(Counter)

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        umi_col = "Combined_UMI" if "Combined_UMI" in (reader.fieldnames or []) else "Combined UMI"
        bc_col = "Barcode"
        for row in reader:
            bc = row[bc_col]
            umi = row[umi_col][: cfg.umi_group_len].upper()

            umi_group = None
            for candidate in bc_umi_groups[bc]:
                if hamming(umi, candidate) is not None and hamming(umi, candidate) <= cfg.umi_group_max_dist:
                    umi_group = candidate
                    break
            if umi_group is None:
                umi_group = umi

            bc_umi_groups[bc][umi_group] += 1

    counts = Counter()
    for bc, umi_groups in bc_umi_groups.items():
        for umi_group, cnt in umi_groups.items():
            counts[f"{bc}_{umi_group}"] += cnt
    return counts


def _split_key(key: str) -> Tuple[str, str]:
    bc, umi = key.rsplit("_", 1)
    return bc, umi


def cluster_barcodes(counts: Counter, cfg: Config) -> Tuple[Counter, Dict[str, str]]:
    """Greedy clustering of barcode variants by Levenshtein distance.

    Returns (rep_bc_counts, bc_to_rep_bc) where rep_bc_counts maps a
    representative ("true") barcode to its aggregated read count, and
    bc_to_rep_bc maps every observed barcode to the representative it was
    collapsed into.
    """
    keys = sorted(counts.keys(), key=lambda k: counts[k], reverse=True)
    key_to_rep: Dict[str, str] = {}

    for key in keys:
        if key in key_to_rep:
            continue
        rep = key
        key_to_rep[key] = rep
        rep_count = counts[rep]
        rep_bc, _ = _split_key(rep)

        for key2 in keys:
            if key2 in key_to_rep:
                continue
            if counts[key2] >= cfg.cluster_protect_count:
                continue
            bc2, _ = _split_key(key2)
            if Levenshtein.distance(rep_bc, bc2) <= cfg.cluster_max_dist:
                if counts[key2] <= cfg.cluster_merge_ratio * rep_count:
                    key_to_rep[key2] = rep

    rep_bc_counts = Counter()
    bc_to_rep_bc: Dict[str, str] = {}
    for key, rep_key in key_to_rep.items():
        bc, _ = _split_key(key)
        rep_bc, _ = _split_key(rep_key)
        rep_bc_counts[rep_bc] += counts[key]
        bc_to_rep_bc[bc] = rep_bc

    return rep_bc_counts, bc_to_rep_bc


def cluster_sample(tsv_path: str, sample_name: str, cfg: Config, out_csv: str) -> Counter:
    counts = load_barcode_umi_counts(tsv_path, cfg)
    rep_bc_counts, _bc_mapping = cluster_barcodes(counts, cfg)

    total = sum(rep_bc_counts.values())
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "corrected_barcode", "total_count", "freq_within_sample"])
        for rep_bc, n in sorted(rep_bc_counts.items(), key=lambda x: -x[1]):
            freq = n / total if total else 0.0
            writer.writerow([sample_name, rep_bc, n, f"{freq:.6f}"])

    return rep_bc_counts
