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
each sample's barcodes can be filtered by a RELATIVE-ABUNDANCE value
threshold (`min_rel_abundance`): drop every barcode whose within-sample
frequency is below `fraction * reference_stat`, where reference_stat is the
sample's mean or median barcode frequency. This is a VALUE cutoff, not a
fixed-fraction cutoff -- it drops however many barcodes fall below the
threshold, whether that is 5%% or 90%% of the barcodes in the file. It
removes the long low-abundance tail (likely sequencing error /
contamination) on a per-sample basis, so a barcode that is abundant in one
sample but only a low-abundance straggler in another is not spuriously
reported as "shared".

Example: with fraction=0.05 and relative_to="median", a barcode is dropped
if its within-sample frequency is below 5%% of the sample's median barcode
frequency.
"""
from __future__ import annotations

import csv
import glob
import os
import statistics
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def _relative_abundance_cutoff(
    counts: Dict[str, int],
    fraction: float,
    relative_to: str = "median",
) -> Tuple[List[str], Optional[float], Optional[float]]:
    """Relative-abundance value threshold, on within-sample frequency.

    Drops barcodes whose within-sample frequency is strictly below
    `fraction * reference_freq`, where reference_freq is the mean or median
    of the sample's per-barcode within-sample frequencies.

    Returns (barcodes_to_drop, cutoff_freq, reference_freq). All three are
    None-safe: if the filter is disabled or cannot be computed, returns
    ([], None, None). cutoff_freq / reference_freq are frequencies (fractions
    of the sample total), used for logging and to draw the plot cutoff line.
    """
    n = len(counts)
    if n <= 1 or fraction <= 0:
        return [], None, None

    total = sum(counts.values())
    if total <= 0:
        return [], None, None

    freqs = {bc: c / total for bc, c in counts.items()}
    values = list(freqs.values())
    if relative_to == "mean":
        reference_freq = statistics.fmean(values)
    elif relative_to == "median":
        reference_freq = statistics.median(values)
    else:
        raise ValueError(f"relative_to must be 'mean' or 'median', got {relative_to!r}")

    cutoff_freq = fraction * reference_freq
    drop = [bc for bc, f in freqs.items() if f < cutoff_freq]
    # deterministic ordering of the dropped list (lowest freq first)
    drop.sort(key=lambda bc: (freqs[bc], bc))
    return drop, cutoff_freq, reference_freq


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
    min_rel_abundance: float = 0.0,
    relative_to: str = "median",
    dist_out: Optional[Dict[str, dict]] = None,
) -> Dict[str, Dict[str, int]]:
    """Returns {barcode: {sample: count}} across all *.csv files in csv_dir.

    min_count:         absolute minimum per-sample count for a barcode
                       (applied first, per barcode).
    min_rel_abundance: if > 0, for each sample independently, drop every
                       barcode whose within-sample frequency is below
                       min_rel_abundance * (mean|median barcode frequency).
                       This is a VALUE threshold: it drops however many
                       barcodes fall below it (possibly a large fraction).
                       Applied after min_count, before pooling across samples.
    relative_to:       "mean" or "median" -- the reference statistic the
                       fraction is applied to. Default "median" (robust to a
                       few dominant barcodes).
    dist_out:          if a dict is provided, it is populated with per-sample
                       distribution info for plotting:
                       {sample: {"kept_counts": [...], "dropped_counts": [...],
                                 "cutoff_freq": float|None,
                                 "reference_freq": float|None,
                                 "relative_to": str,
                                 "sample_total": int, "n_dropped": int,
                                 "n_total": int}}
    """
    barcode_sample_counts: Dict[str, Dict[str, int]] = defaultdict(dict)
    import sys

    for path in sorted(glob.glob(os.path.join(csv_dir, "*.csv"))):
        sample_name, counts = _read_sample_counts(path, min_count)
        if not counts:
            continue

        n_total = len(counts)
        sample_total = sum(counts.values())
        drop_bcs, cutoff_freq, reference_freq = _relative_abundance_cutoff(
            counts, min_rel_abundance, relative_to
        )
        drop_set = set(drop_bcs)
        dropped_counts = [counts[bc] for bc in drop_bcs]
        kept = {bc: c for bc, c in counts.items() if bc not in drop_set}

        if min_rel_abundance and min_rel_abundance > 0:
            if cutoff_freq is not None:
                pct_dropped = 100.0 * len(drop_bcs) / n_total if n_total else 0.0
                print(
                    f"[map] {sample_name}: dropped {len(drop_bcs)}/{n_total} "
                    f"({pct_dropped:.1f}%) barcodes below "
                    f"{min_rel_abundance:g}x{relative_to} "
                    f"({relative_to} freq={reference_freq:.3e}, "
                    f"cutoff freq={cutoff_freq:.3e})",
                    file=sys.stderr,
                )
            else:
                print(
                    f"[map] {sample_name}: relative-abundance filter not applied "
                    f"(only {n_total} barcode(s))",
                    file=sys.stderr,
                )

        if dist_out is not None:
            dist_out[sample_name] = {
                "kept_counts": sorted(kept.values()),
                "dropped_counts": sorted(dropped_counts),
                "cutoff_freq": cutoff_freq,
                "reference_freq": reference_freq,
                "relative_to": relative_to,
                "sample_total": sample_total,
                "n_dropped": len(drop_bcs),
                "n_total": n_total,
            }

        for bc, count in kept.items():
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


def plot_count_distributions(
    dist: Dict[str, dict],
    outdir: str,
    min_rel_abundance: float,
) -> List[str]:
    """Per-sample histogram of log10(within-sample frequency) with the
    relative-abundance cutoff marked.

    Kept vs dropped barcodes are drawn as stacked histograms so you can see
    exactly which part of the tail the cutoff removed. Returns the list of
    written PNG paths. Silently no-ops (with a stderr note) if matplotlib is
    not installed, so the core pipeline never hard-depends on it.
    """
    import sys
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print(
            "[map] matplotlib not available; skipping distribution plots. "
            "Install it (conda install -c conda-forge matplotlib) to enable.",
            file=sys.stderr,
        )
        return []

    import math as _math

    plot_dir = os.path.join(outdir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    written: List[str] = []

    for sample, info in sorted(dist.items()):
        total = info["sample_total"] or 1
        # Convert counts -> within-sample frequency, then log10.
        kept_freq = [c / total for c in info["kept_counts"]]
        dropped_freq = [c / total for c in info["dropped_counts"]]
        kept_log = [_math.log10(f) for f in kept_freq if f > 0]
        dropped_log = [_math.log10(f) for f in dropped_freq if f > 0]

        all_log = kept_log + dropped_log
        if not all_log:
            continue

        lo, hi = min(all_log), max(all_log)
        if lo == hi:
            lo -= 0.5
            hi += 0.5
        n_bins = max(10, min(60, int(_math.sqrt(info["n_total"]) * 2)))
        bins = [lo + (hi - lo) * i / n_bins for i in range(n_bins + 1)]

        fig, ax = plt.subplots(figsize=(8, 5))
        if dropped_log:
            ax.hist(
                [kept_log, dropped_log], bins=bins, stacked=True,
                color=["#4C78A8", "#E45756"],
                label=[f"kept (n={len(kept_log)})", f"dropped (n={len(dropped_log)})"],
            )
        else:
            ax.hist(kept_log, bins=bins, color="#4C78A8", label=f"kept (n={len(kept_log)})")

        # Cutoff line at fraction * reference_freq.
        cutoff_freq = info.get("cutoff_freq")
        reference_freq = info.get("reference_freq")
        relative_to = info.get("relative_to", "median")
        if cutoff_freq is not None and cutoff_freq > 0:
            x = _math.log10(cutoff_freq)
            ax.axvline(
                x, color="black", linestyle="--", linewidth=1.5,
                label=(
                    f"cutoff = {min_rel_abundance:g}×{relative_to}\n"
                    f"(freq < {cutoff_freq:.2e})"
                ),
            )
        # Reference-stat line (mean/median) for context.
        if reference_freq is not None and reference_freq > 0:
            ax.axvline(
                _math.log10(reference_freq), color="gray", linestyle=":", linewidth=1.2,
                label=f"{relative_to} freq ({reference_freq:.2e})",
            )

        pct_dropped = 100.0 * info["n_dropped"] / info["n_total"] if info["n_total"] else 0.0
        ax.set_xlabel("log10(within-sample frequency)")
        ax.set_ylabel("number of barcodes")
        ax.set_title(
            f"{sample} — barcode frequency distribution\n"
            f"{info['n_total']} barcodes, dropped {info['n_dropped']} "
            f"({pct_dropped:.1f}%) below {min_rel_abundance:g}×{relative_to}"
        )
        ax.legend(fontsize=8)
        fig.tight_layout()

        out_png = os.path.join(plot_dir, f"{sample}.count_distribution.png")
        fig.savefig(out_png, dpi=120)
        plt.close(fig)
        written.append(out_png)

    return written
