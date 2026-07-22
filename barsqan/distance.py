"""
Distance helpers with a robust backend.

rapidfuzz ships compiled C extensions for most CPython versions, but on very
new interpreters (e.g. Python 3.14 before wheels are published) it silently
falls back to a pure-Python implementation that is ~100x slower. When the
clustering step calls it millions of times, that fallback makes the pipeline
appear to hang.

To stay fast everywhere we:
  * use rapidfuzz's compiled backend when it is actually the C extension;
  * otherwise use our own bounded (banded) Levenshtein with an early-exit
    score cutoff, which is more than fast enough for short barcode strings
    and never depends on a compiled extension being present.

Both paths honor a `max_dist` cutoff so we can stop as soon as the distance
is known to exceed the clustering threshold.
"""
from __future__ import annotations

from typing import Optional


def hamming(s1: str, s2: str) -> Optional[int]:
    """Hamming distance for equal-length strings; None if lengths differ."""
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


# ---- Decide, once at import time, whether rapidfuzz's C backend is usable ----
_USE_RAPIDFUZZ_CPP = False
try:  # pragma: no cover - depends on install environment
    from rapidfuzz.distance import Levenshtein as _RF_Lev

    _mod = getattr(_RF_Lev.distance, "__module__", "")
    # The compiled backend lives in rapidfuzz.distance.metrics_cpp* ; the slow
    # pure-Python fallback lives in ...Levenshtein_py. Only trust the former.
    _USE_RAPIDFUZZ_CPP = "cpp" in _mod
except Exception:
    _USE_RAPIDFUZZ_CPP = False


def _bounded_levenshtein(a: str, b: str, max_dist: int) -> Optional[int]:
    """Levenshtein distance if <= max_dist, else None. Banded DP, O(n*k)."""
    la, lb = len(a), len(b)
    if abs(la - lb) > max_dist:
        return None
    if la == 0:
        return lb if lb <= max_dist else None
    if lb == 0:
        return la if la <= max_dist else None

    # Ensure a is the shorter string to keep the row small.
    if la > lb:
        a, b = b, a
        la, lb = lb, la

    INF = max_dist + 1
    prev = list(range(la + 1))
    for j in range(1, lb + 1):
        cur = [j] + [0] * la
        row_min = j
        bj = b[j - 1]
        # Band: only i within max_dist of j can matter.
        lo = max(1, j - max_dist)
        hi = min(la, j + max_dist)
        if lo > 1:
            cur[lo - 1] = INF
        for i in range(lo, hi + 1):
            cost = 0 if a[i - 1] == bj else 1
            v = min(
                prev[i] + 1,        # deletion
                cur[i - 1] + 1,     # insertion
                prev[i - 1] + cost, # substitution
            )
            cur[i] = v
            if v < row_min:
                row_min = v
        if hi < la:
            cur[hi + 1] = INF
        if row_min > max_dist:
            return None
        prev = cur

    d = prev[la]
    return d if d <= max_dist else None


def levenshtein_within(a: str, b: str, max_dist: int) -> Optional[int]:
    """Return Levenshtein distance if it is <= max_dist, else None.

    Uses rapidfuzz's compiled backend with score_cutoff when available,
    otherwise a self-contained bounded DP. The None-vs-value contract is
    identical across backends so callers need not care which is active.
    """
    if _USE_RAPIDFUZZ_CPP:
        # rapidfuzz returns max_dist+1 (or the configured score_cutoff+1
        # sentinel) when the distance exceeds the cutoff.
        d = _RF_Lev.distance(a, b, score_cutoff=max_dist)
        return d if d <= max_dist else None
    return _bounded_levenshtein(a, b, max_dist)


def backend_name() -> str:
    return "rapidfuzz-cpp" if _USE_RAPIDFUZZ_CPP else "pure-python-bounded"
