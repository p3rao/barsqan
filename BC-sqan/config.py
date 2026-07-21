"""
Central configuration for the barcode toolkit.

All amplicon-design parameters (primers, index/UMI lengths, barcode motif,
quality thresholds, overlap requirements, clustering tolerances) live here in
one place so extraction / filtering / clustering all stay consistent with a
single source of truth, instead of being hardcoded separately in each script.

Load a YAML file with `Config.from_yaml(path)` to override any subset of the
defaults below. Any keys not present in the YAML keep their default value.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional

try:
    import yaml  # PyYAML
except ImportError:  # pragma: no cover
    yaml = None


@dataclass
class Config:
    # ---------------- Read layout ----------------
    # Read 1: INDEX1 + UMI1 + INNER_P1_PRIMER + BARCODE + DOWNSTREAM1
    # Read 2 (paired-end only): INDEX2 + UMI2 + INNER_P2_PRIMER + DOWNSTREAM2
    index1_len: int = 6
    umi1_len: int = 8
    index2_len: int = 6
    umi2_len: int = 8

    inner_p1_primer: str = "GCCACTCGAGCACCTAGGAGG"
    inner_p2_primer: str = "AATCAGTCCAGTTATGCTGTGAA"

    # Max edit distance (substitutions/indels) tolerated when fuzzy-matching
    # each primer (uses the `regex` module's {e<=k} fuzzy matching).
    primer1_max_edits: int = 3
    primer2_max_edits: int = 3

    # ---------------- Barcode motif ----------------
    # Structured barcode motif, e.g. "NNNNNacNNNNNgaNNNNNtcNNNNN".
    # Uppercase 'N' = any base (wildcard). Lowercase letters = fixed bases.
    # This single string drives both the strict extraction regex AND the
    # tolerant motif re-check used during filtering, so the two steps can
    # never drift out of sync with each other.
    barcode_motif: str = "NNNNNacNNNNNgaNNNNNtcNNNNN"

    # ---------------- Paired-end overlap validation ----------------
    # Minimum length of the downstream constant region required after the
    # barcode (R1) / after primer2 (R2) for the read to be considered
    # structurally complete.
    downstream_min_len: int = 20
    # Sliding-window size used to find the best-aligned overlap between the
    # R1 downstream sequence and the reverse complement of the R2 downstream
    # sequence (handles the fact that R1/R2 downstream lengths need not be
    # positionally aligned).
    overlap_k: int = 20
    # Max Hamming mismatches allowed within the best k-mer alignment.
    overlap_max_mismatch: int = 4

    # ---------------- Quality filtering ----------------
    # Minimum Phred quality (Phred+33) required across the *structured*
    # prefix of each read (index+UMI+primer[+barcode for R1]). We deliberately
    # do NOT require high quality across the long constant tail - that
    # region is only used for structural/overlap validation, and requiring
    # Q30+ there discards the majority of real reads for no benefit.
    qmin: int = 25

    # ---------------- Internal-index filtering ----------------
    # Length of each index "chunk" checked independently against the
    # expected index (e.g. 6bp index1, optionally 6bp index2 for PE => 12bp
    # combined). Each chunk independently tolerates up to index_max_mismatch.
    index_chunk_len: int = 6
    index_max_mismatch: int = 1

    # Tolerance used when re-checking the barcode motif during filtering
    # (independent, tolerant re-verification of the same barcode_motif).
    motif_max_mismatch: int = 1

    # ---------------- UMI collapsing / barcode clustering ----------------
    # Number of bases (from the start of the combined UMI) used to group
    # near-identical UMIs together (Hamming distance <= umi_group_max_dist).
    umi_group_len: int = 8
    umi_group_max_dist: int = 1

    # Barcode clustering (collapsing sequencing-error variants of the same
    # true barcode) via Levenshtein distance.
    cluster_max_dist: int = 2
    # A low-count barcode variant is merged into a higher-count "parent"
    # barcode only if its count is <= cluster_merge_ratio * parent_count.
    cluster_merge_ratio: float = 0.10
    # Barcodes with count >= cluster_protect_count are never merged into
    # another barcode (they are always treated as their own true barcode).
    cluster_protect_count: int = 25

    # ---------------- Misc ----------------
    # How many of the most frequent regex-non-matching R1 sequences to save
    # per sample (helps spot off-target / contaminating amplicons).
    top_unmatched_n: int = 50

    @classmethod
    def from_yaml(cls, path: Optional[str]) -> "Config":
        cfg = cls()
        if not path:
            return cfg
        if yaml is None:
            raise RuntimeError("PyYAML is required to load a config file (pip install pyyaml)")
        with open(path) as f:
            overrides = yaml.safe_load(f) or {}
        valid_fields = {f.name for f in dataclasses.fields(cls)}
        unknown = set(overrides) - valid_fields
        if unknown:
            raise ValueError(f"Unknown config key(s) in {path}: {sorted(unknown)}")
        return dataclasses.replace(cfg, **overrides)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)
