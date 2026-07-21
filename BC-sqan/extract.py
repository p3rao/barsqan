"""
Step 1: extract barcodes (+ index/UMI) from raw FASTQ reads.

Read layout enforced (see config.py for all tunables):
  R1 (always):            INDEX1 + UMI1 + INNER_P1_PRIMER(fuzzy) + BARCODE(motif) + DOWNSTREAM1
  R2 (paired-end only):    INDEX2 + UMI2 + INNER_P2_PRIMER(fuzzy) + DOWNSTREAM2

For paired-end data, DOWNSTREAM1 (R1) and the reverse complement of
DOWNSTREAM2 (R2) are expected to overlap in the shared constant region
downstream of the barcode; this is used as a structural sanity check
(best-aligned sliding-window Hamming distance) before a read pair is
accepted.

For single-end data, only the R1 structure above is required (no overlap
check is possible/needed).

Any read (or pair) that fails the quality filter or the structural regex
is counted, and for regex failures the raw R1 sequence is tallied so you
can see what fraction of "junk" reads come from a handful of recurring
off-target/contaminating amplicons.
"""
from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Optional

import regex  # pip install regex

from .config import Config
from .fastq_io import SampleFiles, read_fastq
from .motif import motif_to_regex, motif_length


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming_distance(s1: str, s2: str) -> Optional[int]:
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def best_overlap_distance(s1: str, s2: str, k: int) -> Optional[int]:
    """Minimal Hamming distance between any k-mer of s1 and any k-mer of s2."""
    s1u, s2u = s1.upper(), s2.upper()
    if len(s1u) < k or len(s2u) < k:
        return None
    best = None
    for i in range(len(s1u) - k + 1):
        w1 = s1u[i:i + k]
        for j in range(len(s2u) - k + 1):
            w2 = s2u[j:j + k]
            d = hamming_distance(w1, w2)
            if best is None or d < best:
                best = d
                if best == 0:
                    return 0
    return best


def region_all_q(qual_str: str, start: int, length: int, qmin: int) -> bool:
    end = start + length
    if len(qual_str) < end:
        return False
    threshold = qmin + 33  # Phred+33
    return all(ord(c) >= threshold for c in qual_str[start:end])


@dataclass
class ExtractStats:
    total_reads: int = 0
    quality_fail: int = 0
    regex_fail: int = 0
    overlap_fail: int = 0
    passed: int = 0

    def pass_fraction(self) -> float:
        return self.passed / self.total_reads if self.total_reads else 0.0

    def to_lines(self) -> str:
        return (
            f"# Total_reads\t{self.total_reads}\n"
            f"# Passed\t{self.passed}\n"
            f"# Quality_fail\t{self.quality_fail}\n"
            f"# Regex_fail\t{self.regex_fail}\n"
            f"# Overlap_fail\t{self.overlap_fail}\n"
            f"# Pass_fraction\t{self.pass_fraction():.4f}\n"
        )


class AmpliconExtractor:
    """Builds the R1/R2 regexes once from a Config and reuses them for every read."""

    def __init__(self, cfg: Config):
        self.cfg = cfg
        bc_group = motif_to_regex(cfg.barcode_motif)
        bc_len = motif_length(cfg.barcode_motif)

        idx1 = f"(.{{{cfg.index1_len}}})"
        umi1 = f"(.{{{cfg.umi1_len}}})"
        idx2 = f"(.{{{cfg.index2_len}}})"
        umi2 = f"(.{{{cfg.umi2_len}}})"

        # R1: idx1(g1) umi1(g2) primer1(fuzzy,noncapture) barcode(g3) downstream(g4)
        r1_pattern = (
            r"^" + idx1 + umi1
            + r"(?:%s){e<=%d}" % (cfg.inner_p1_primer, cfg.primer1_max_edits)
            + bc_group
            + r"(" + r"[ACGTN]{%d,}" % cfg.downstream_min_len + r")"
            + r".*$"
        )
        self.read1_regex = regex.compile(r1_pattern, regex.IGNORECASE)

        # R2 (PE only): idx2(g1) umi2(g2) primer2(fuzzy,noncapture) downstream(g3)
        r2_pattern = (
            r"^" + idx2 + umi2
            + r"(?:%s){e<=%d}" % (cfg.inner_p2_primer, cfg.primer2_max_edits)
            + r"(" + r"[ACGTN]{%d,}" % cfg.downstream_min_len + r")"
            + r".*$"
        )
        self.read2_regex = regex.compile(r2_pattern, regex.IGNORECASE)

        self.bc_len = bc_len
        self.r1_qc_len = cfg.index1_len + cfg.umi1_len + len(cfg.inner_p1_primer) + bc_len
        self.r2_qc_len = cfg.index2_len + cfg.umi2_len + len(cfg.inner_p2_primer)

    # -------------------- paired-end --------------------
    def process_pair(self, r1_seq, r1_qual, r2_seq, r2_qual, stats: ExtractStats, unmatched: Counter):
        cfg = self.cfg
        stats.total_reads += 1

        if not region_all_q(r1_qual, 0, self.r1_qc_len, cfg.qmin) or \
           not region_all_q(r2_qual, 0, self.r2_qc_len, cfg.qmin):
            stats.quality_fail += 1
            return None

        m1 = self.read1_regex.match(r1_seq)
        m2 = self.read2_regex.match(r2_seq)
        if not (m1 and m2):
            stats.regex_fail += 1
            unmatched[r1_seq[:60]] += 1
            return None

        idx1, umi1, barcode, r1_down = m1.group(1), m1.group(2), m1.group(3), m1.group(4)
        idx2, umi2, r2_down = m2.group(1), m2.group(2), m2.group(3)

        r2_down_rc = rev_comp(r2_down)
        best_dist = best_overlap_distance(r1_down, r2_down_rc, cfg.overlap_k)
        if best_dist is None or best_dist > cfg.overlap_max_mismatch:
            stats.overlap_fail += 1
            return None

        stats.passed += 1
        return {
            "Combined_UMI": umi1 + umi2,
            "Combined_Index": idx1 + idx2,
            "Barcode": barcode,
            "Index1": idx1,
            "Index2": idx2,
            "UMI1": umi1,
            "UMI2": umi2,
            "BestOverlapDist": best_dist,
        }

    # -------------------- single-end --------------------
    def process_single(self, r1_seq, r1_qual, stats: ExtractStats, unmatched: Counter):
        cfg = self.cfg
        stats.total_reads += 1

        if not region_all_q(r1_qual, 0, self.r1_qc_len, cfg.qmin):
            stats.quality_fail += 1
            return None

        m1 = self.read1_regex.match(r1_seq)
        if not m1:
            stats.regex_fail += 1
            unmatched[r1_seq[:60]] += 1
            return None

        idx1, umi1, barcode, _r1_down = m1.group(1), m1.group(2), m1.group(3), m1.group(4)
        stats.passed += 1
        return {
            "Combined_UMI": umi1,
            "Combined_Index": idx1,
            "Barcode": barcode,
            "Index1": idx1,
            "Index2": "",
            "UMI1": umi1,
            "UMI2": "",
            "BestOverlapDist": "",
        }


PE_COLUMNS = ["Combined_UMI", "Combined_Index", "Barcode", "Index1", "Index2", "UMI1", "UMI2", "BestOverlapDist"]
SE_COLUMNS = ["Combined_UMI", "Combined_Index", "Barcode", "Index1", "UMI1"]


def extract_sample(sf: SampleFiles, cfg: Config, outdir: str) -> ExtractStats:
    """Run extraction for one sample and write:
      {outdir}/{sample}.parsed.tsv        - passing reads
      {outdir}/{sample}.stats.tsv         - pass/fail counters
      {outdir}/{sample}.unmatched.tsv     - top N most frequent regex-non-matching R1 sequences
    """
    os.makedirs(outdir, exist_ok=True)
    extractor = AmpliconExtractor(cfg)
    stats = ExtractStats()
    unmatched: Counter = Counter()

    parsed_path = os.path.join(outdir, f"{sf.sample_name}.parsed.tsv")
    columns = PE_COLUMNS if sf.mode == "PE" else SE_COLUMNS

    with open(parsed_path, "w") as out:
        out.write("\t".join(columns) + "\n")

        if sf.mode == "PE":
            r1_iter = read_fastq(sf.r1)
            r2_iter = read_fastq(sf.r2)
            for (h1, s1, q1), (h2, s2, q2) in zip(r1_iter, r2_iter):
                row = extractor.process_pair(s1, q1, s2, q2, stats, unmatched)
                if row is not None:
                    out.write("\t".join(str(row[c]) for c in columns) + "\n")
        else:
            for (h1, s1, q1) in read_fastq(sf.r1):
                row = extractor.process_single(s1, q1, stats, unmatched)
                if row is not None:
                    out.write("\t".join(str(row[c]) for c in columns) + "\n")

    stats_path = os.path.join(outdir, f"{sf.sample_name}.stats.tsv")
    with open(stats_path, "w") as f:
        f.write(stats.to_lines())

    unmatched_path = os.path.join(outdir, f"{sf.sample_name}.unmatched.tsv")
    with open(unmatched_path, "w") as f:
        f.write("sequence_prefix\tcount\n")
        for seq, cnt in unmatched.most_common(cfg.top_unmatched_n):
            f.write(f"{seq}\t{cnt}\n")

    return stats
