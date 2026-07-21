"""
Shared barcode-motif logic.

A single motif string like "NNNNNacNNNNNgaNNNNNtcNNNNN" drives both:
  1. the STRICT extraction regex (used while parsing FASTQ reads), and
  2. the TOLERANT re-check used during filtering (allows a few mismatches
     at the fixed positions, in case of sequencing error).

Keeping both derived from one motif string means they can never silently
drift apart the way they did across the original standalone scripts.
"""
from __future__ import annotations

from typing import List, Tuple


def motif_length(motif: str) -> int:
    return len(motif)


def motif_to_regex(motif: str) -> str:
    """Build a capturing regex group that strictly matches `motif`.

    Uppercase 'N' -> any base ([ACGTN]). Any other (lowercase) character is
    treated as a fixed literal base at that position.
    """
    parts = []
    for ch in motif:
        if ch == "N":
            parts.append("[ACGTN]")
        elif ch.upper() in "ACGT":
            parts.append(ch.upper())
        else:
            raise ValueError(f"Unsupported character '{ch}' in barcode motif '{motif}'")
    return "(" + "".join(parts) + ")"


def motif_const_positions(motif: str) -> List[Tuple[int, str]]:
    """Return [(position, expected_base_uppercase), ...] for fixed (non-N) positions."""
    out = []
    for i, ch in enumerate(motif):
        if ch != "N":
            out.append((i, ch.upper()))
    return out


def barcode_matches_motif(barcode: str, motif: str, max_mismatch: int = 1) -> bool:
    """Tolerant re-check: does `barcode` contain a window matching `motif`
    (allowing up to max_mismatch mismatches at the fixed positions)?

    Slides a window of len(motif) across `barcode` in case extraction
    captured a slightly longer/shorter downstream-inclusive fragment.
    """
    b = barcode.upper()
    m_len = len(motif)
    const_positions = motif_const_positions(motif)
    if len(b) < m_len:
        return False
    for start in range(0, len(b) - m_len + 1):
        window = b[start:start + m_len]
        mismatches = 0
        for pos, expected_base in const_positions:
            if window[pos] != expected_base:
                mismatches += 1
                if mismatches > max_mismatch:
                    break
        if mismatches <= max_mismatch:
            return True
    return False
