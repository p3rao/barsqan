"""
FASTQ reading and sample-name -> file discovery.

Handles:
- plain and gzip-compressed FASTQ (.fastq/.fq, .fastq.gz/.fq.gz)
- matching a list of user-supplied sample names against files in a directory
- auto-detecting single-end (1 file per sample) vs paired-end (R1+R2) layout
"""
from __future__ import annotations

import gzip
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Tuple

FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Common ways R1/R2 are encoded in Illumina-style filenames.
_R1_PATTERNS = [r"_R1[_.]", r"_R1$", r"_1\.fastq", r"_1\.fq"]
_R2_PATTERNS = [r"_R2[_.]", r"_R2$", r"_2\.fastq", r"_2\.fq"]


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fastq(path: str) -> Iterator[Tuple[str, str, str]]:
    """Yield (header, sequence, quality) tuples from a FASTQ file."""
    with open_maybe_gzip(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip("\n")
            plus = f.readline()
            qual = f.readline().rstrip("\n")
            if not plus:
                raise ValueError(f"Truncated FASTQ record in {path}")
            yield header.rstrip("\n"), seq, qual


def _is_r1(filename: str) -> bool:
    return any(re.search(p, filename) for p in _R1_PATTERNS)


def _is_r2(filename: str) -> bool:
    return any(re.search(p, filename) for p in _R2_PATTERNS)


def _list_fastq_files(fastq_dir: str) -> List[str]:
    files = []
    for fn in os.listdir(fastq_dir):
        if fn.lower().endswith(FASTQ_EXTS):
            files.append(fn)
    return sorted(files)


def _sample_matches_filename(sample_name: str, filename: str) -> bool:
    """True if sample_name appears in filename at a token boundary.

    Avoids partial-name collisions like "Col-1" matching "Col-10" by
    requiring the characters immediately before/after the match (if any)
    to be non-alphanumeric.
    """
    pattern = r"(?<![A-Za-z0-9])" + re.escape(sample_name) + r"(?![A-Za-z0-9])"
    return re.search(pattern, filename) is not None


@dataclass
class SampleFiles:
    sample_name: str
    mode: str  # "SE" or "PE"
    r1: str
    r2: Optional[str] = None


def find_sample_files(
    sample_names: List[str],
    fastq_dir: str,
    mode: str = "auto",
) -> Dict[str, SampleFiles]:
    """Match each sample name to its FASTQ file(s) in fastq_dir.

    mode: "auto" (detect per sample from file count), "se", or "pe".

    Returns dict {sample_name: SampleFiles}. Raises ValueError with a
    descriptive message for samples with no match or an ambiguous match
    that can't be resolved.
    """
    all_files = _list_fastq_files(fastq_dir)
    # Longest sample names first, so a longer/more-specific name is matched
    # before a shorter name that happens to be a prefix of it.
    ordered_samples = sorted(sample_names, key=len, reverse=True)

    result: Dict[str, SampleFiles] = {}
    used_files = set()

    for sample in ordered_samples:
        candidates = [
            fn for fn in all_files
            if fn not in used_files and _sample_matches_filename(sample, fn)
        ]
        if not candidates:
            raise ValueError(
                f"No FASTQ files found for sample '{sample}' in {fastq_dir}. "
                f"Checked {len(all_files)} files."
            )

        r1_candidates = [fn for fn in candidates if _is_r1(fn)]
        r2_candidates = [fn for fn in candidates if _is_r2(fn)]
        other_candidates = [fn for fn in candidates if fn not in r1_candidates and fn not in r2_candidates]

        forced_mode = None if mode == "auto" else mode.upper()

        if r1_candidates and r2_candidates:
            if len(r1_candidates) > 1 or len(r2_candidates) > 1:
                raise ValueError(
                    f"Ambiguous match for sample '{sample}': "
                    f"R1 candidates={r1_candidates} R2 candidates={r2_candidates}"
                )
            if forced_mode == "SE":
                raise ValueError(
                    f"Sample '{sample}' has both R1 and R2 files but mode=se was forced."
                )
            sf = SampleFiles(
                sample_name=sample, mode="PE",
                r1=os.path.join(fastq_dir, r1_candidates[0]),
                r2=os.path.join(fastq_dir, r2_candidates[0]),
            )
            used_files.update(r1_candidates + r2_candidates)
        elif len(other_candidates) == 1 and not r1_candidates and not r2_candidates:
            if forced_mode == "PE":
                raise ValueError(
                    f"Sample '{sample}' matched only one file but mode=pe was forced: "
                    f"{other_candidates[0]}"
                )
            sf = SampleFiles(
                sample_name=sample, mode="SE",
                r1=os.path.join(fastq_dir, other_candidates[0]),
                r2=None,
            )
            used_files.update(other_candidates)
        elif len(candidates) == 1:
            # A single file that happens to look like an R1 (no R2 counterpart) - SE.
            sf = SampleFiles(
                sample_name=sample, mode="SE",
                r1=os.path.join(fastq_dir, candidates[0]),
                r2=None,
            )
            used_files.update(candidates)
        else:
            raise ValueError(
                f"Could not unambiguously resolve FASTQ file(s) for sample '{sample}'. "
                f"Candidates found: {candidates}"
            )
        result[sample] = sf

    return result
