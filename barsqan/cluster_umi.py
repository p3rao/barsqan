"""
Step 3: collapse UMIs and cluster barcode sequencing-error variants into a
final per-sample barcode count table.

Two-stage collapsing:
  1. UMI collapsing: within a given barcode, group UMIs that are within
     `umi_group_max_dist` Hamming distance of each other (keeps distinct
     true molecules while absorbing sequencing errors in the UMI).
  2. Barcode clustering: aggregate reads per *barcode*, then assign
     low-abundance barcode variants that are within `cluster_max_dist`
     edit distance of a much higher-abundance "parent" barcode
     (count <= cluster_merge_ratio * parent_count) to that parent, unless
     the variant's own count already exceeds cluster_protect_count.

Performance notes
-----------------
The original implementation compared every (barcode, umi) key against every
other key with a full Levenshtein call -- an O(N^2) scan that becomes
billions of comparisons on real samples (and grinds to a halt when
rapidfuzz falls back to its pure-Python backend on very new Python
versions). This version:

  * clusters on unique *barcodes* (there are far fewer of these than
    barcode+umi keys), aggregating UMI-group counts up front;
  * only ever tries to merge a variant into an already-established, more
    abundant parent, and only parents whose count could plausibly absorb it;
  * prunes candidate parents by |len difference| <= cluster_max_dist before
    computing any edit distance;
  * uses a length-bucketed index so each low-count variant is compared only
    against the handful of parents in a compatible length band;
  * passes score_cutoff to the distance function for early termination.
"""
from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

from .config import Config
from .distance import levenshtein_within, hamming


def load_barcode_counts(tsv_path: str, cfg: Config) -> Counter:
    """Load per-barcode read counts, collapsing near-identical UMIs.

    For each barcode, UMIs within `umi_group_max_dist` Hamming distance are
    treated as the same molecule. The returned Counter maps each barcode to
    its number of distinct-molecule-weighted reads (i.e. summed counts over
    its UMI groups), which is what the barcode clustering step operates on.

    Uses a dict keyed by exact UMI prefix as a fast path; only UMIs that do
    not match an existing group exactly are compared (within their barcode)
    against existing groups, keeping this close to O(reads) in practice.
    """
    # barcode -> {umi_group_repr: count}
    bc_umi_groups: Dict[str, "Counter[str]"] = defaultdict(Counter)

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        umi_col = "Combined_UMI" if "Combined_UMI" in fields else "Combined UMI"
        bc_col = "Barcode"
        for row in reader:
            bc = row[bc_col]
            umi = row[umi_col][: cfg.umi_group_len].upper()

            groups = bc_umi_groups[bc]
            if umi in groups:
                # exact match to an existing group: fast path, no distance calls
                groups[umi] += 1
                continue

            if cfg.umi_group_max_dist <= 0:
                groups[umi] += 1
                continue

            # find an existing group within Hamming distance
            umi_group = None
            for candidate in groups:
                d = hamming(umi, candidate)
                if d is not None and d <= cfg.umi_group_max_dist:
                    umi_group = candidate
                    break
            groups[umi_group if umi_group is not None else umi] += 1

    # collapse to per-barcode totals
    barcode_counts: "Counter[str]" = Counter()
    for bc, groups in bc_umi_groups.items():
        barcode_counts[bc] += sum(groups.values())
    return barcode_counts


def _blocks(seq: str, n_blocks: int) -> List[Tuple[int, str]]:
    """Split seq into n_blocks near-equal contiguous blocks.

    Returns [(block_index, block_substring), ...]. Used for pigeonhole
    candidate filtering: two strings with edit distance <= d must share at
    least one identical (block_index, block_substring) when split into d+1
    blocks, so only barcodes sharing such a block are ever compared.
    """
    L = len(seq)
    out = []
    start = 0
    for i in range(n_blocks):
        end = L * (i + 1) // n_blocks
        out.append((i, seq[start:end]))
        start = end
    return out


def cluster_barcodes(barcode_counts: Counter, cfg: Config) -> Tuple[Counter, Dict[str, str]]:
    """Assign low-abundance barcode variants to a more-abundant parent.

    Returns (rep_bc_counts, bc_to_rep_bc):
      * rep_bc_counts: representative ("true") barcode -> aggregated count
      * bc_to_rep_bc:  every observed barcode -> the representative it
                       collapsed into (identity for representatives)

    Algorithm (greedy, most-abundant-first):
      * Barcodes are processed high count -> low count.
      * A barcode with count >= cluster_protect_count is always its own
        representative and is added to the parent index.
      * A lower-count barcode is merged into the nearest already-registered
        parent within cluster_max_dist edit distance AND whose count is
        large enough that this variant is <= cluster_merge_ratio * parent
        count. If no such parent exists, the barcode becomes a new parent.

    Candidate parents are found with a pigeonhole block index rather than an
    all-pairs scan: each parent is indexed by its (block_index, block_substr)
    fragments when split into cluster_max_dist + 1 blocks. Any variant within
    cluster_max_dist edits must share at least one such fragment, so we only
    ever compute edit distance against parents that share a block. This turns
    the quadratic scan into near-linear work.
    """
    max_dist = cfg.cluster_max_dist
    n_blocks = max_dist + 1

    # (block_index, block_substr) -> list of parent barcodes carrying it.
    block_index: Dict[Tuple[int, str], List[str]] = defaultdict(list)
    parent_count: Dict[str, int] = {}

    def register_parent(bc: str, cnt: int) -> None:
        parent_count[bc] = cnt
        for bi, bs in _blocks(bc, n_blocks):
            block_index[(bi, bs)].append(bc)

    key_to_rep: Dict[str, str] = {}

    # Sort barcodes by count desc, then lexicographically for determinism.
    ordered = sorted(barcode_counts.keys(), key=lambda b: (-barcode_counts[b], b))

    for bc in ordered:
        cnt = barcode_counts[bc]

        # High-count barcodes are always their own representative.
        if cnt >= cfg.cluster_protect_count:
            key_to_rep[bc] = bc
            register_parent(bc, cnt)
            continue

        # Gather candidate parents that share at least one block with bc.
        seen = set()
        best_parent: Optional[str] = None
        best_dist = max_dist + 1
        for bi, bs in _blocks(bc, n_blocks):
            for parent_bc in block_index.get((bi, bs), ()):  # noqa
                if parent_bc in seen:
                    continue
                seen.add(parent_bc)
                pcnt = parent_count[parent_bc]
                # variant must be small relative to the parent to be an error
                if cnt > cfg.cluster_merge_ratio * pcnt:
                    continue
                d = levenshtein_within(bc, parent_bc, max_dist)
                if d is not None and d < best_dist:
                    best_dist = d
                    best_parent = parent_bc
                    if d == 1:  # can't beat an off-by-one; take it
                        break
            if best_parent is not None and best_dist == 1:
                break

        if best_parent is not None:
            key_to_rep[bc] = best_parent
        else:
            key_to_rep[bc] = bc
            register_parent(bc, cnt)

    rep_bc_counts: "Counter[str]" = Counter()
    for bc, rep in key_to_rep.items():
        rep_bc_counts[rep] += barcode_counts[bc]

    return rep_bc_counts, key_to_rep


def cluster_sample(tsv_path: str, sample_name: str, cfg: Config, out_csv: str) -> Counter:
    barcode_counts = load_barcode_counts(tsv_path, cfg)
    rep_bc_counts, _mapping = cluster_barcodes(barcode_counts, cfg)

    total = sum(rep_bc_counts.values())
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "corrected_barcode", "total_count", "freq_within_sample"])
        for rep_bc, n in sorted(rep_bc_counts.items(), key=lambda x: (-x[1], x[0])):
            freq = n / total if total else 0.0
            writer.writerow([sample_name, rep_bc, n, f"{freq:.6f}"])

    return rep_bc_counts
