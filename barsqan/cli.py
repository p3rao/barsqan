"""
Command-line interface for barsqan.

Subcommands:
  extract   - find sample FASTQ files, extract index/UMI/barcode -> {sample}.parsed.tsv
  filter    - filter parsed barcodes by expected internal index + motif re-check
  cluster   - collapse UMIs and cluster barcode sequencing-error variants
  map       - build a barcode x sample overlap matrix from cluster outputs
  run       - convenience pipeline: extract -> filter (optional) -> cluster -> map (optional)

See README.md for a full walkthrough and `config.example.yaml` for every
tunable parameter.
"""
from __future__ import annotations

import argparse
import os
import sys

from .config import Config
from .fastq_io import find_sample_files
from .extract import extract_sample
from .filter_reads import filter_sample, load_expected_index_map, write_summaries, FilterCounts
from .cluster_umi import cluster_sample
from .map_overlap import load_all_counts, write_overlap_outputs


def _read_sample_names(path: str):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def cmd_extract(args):
    cfg = Config.from_yaml(args.config)
    sample_names = _read_sample_names(args.samples)
    sample_files = find_sample_files(sample_names, args.fastq_dir, mode=args.mode)

    os.makedirs(args.outdir, exist_ok=True)
    for sample, sf in sample_files.items():
        print(f"[extract] {sample} ({sf.mode}): R1={sf.r1}" + (f" R2={sf.r2}" if sf.r2 else ""), file=sys.stderr)
        stats = extract_sample(sf, cfg, args.outdir)
        print(f"[extract]   -> passed {stats.passed}/{stats.total_reads} ({stats.pass_fraction():.2%})", file=sys.stderr)


def cmd_filter(args):
    cfg = Config.from_yaml(args.config)
    expected_map = load_expected_index_map(args.index_map)

    per_sample_counts = {}
    for parsed_tsv in sorted(_glob_parsed(args.parsed_dir)):
        sample = os.path.basename(parsed_tsv).replace(".parsed.tsv", "")
        if sample not in expected_map:
            print(f"[filter] WARNING: no expected index for sample '{sample}', skipping", file=sys.stderr)
            continue
        counts = filter_sample(
            parsed_tsv, sample, expected_map[sample], cfg,
            kept_out_dir=os.path.join(args.outdir, "kept"),
            filtered_out_dir=os.path.join(args.outdir, "filtered"),
        )
        per_sample_counts[sample] = counts
        print(f"[filter] {sample}: kept={counts.kept} index_mismatch={counts.index_mismatch} motif_absent={counts.motif_absent}", file=sys.stderr)

    write_summaries(per_sample_counts, os.path.join(args.outdir, "filter_summary.csv"))


def _glob_parsed(parsed_dir):
    import glob
    return glob.glob(os.path.join(parsed_dir, "*.parsed.tsv"))


def cmd_cluster(args):
    cfg = Config.from_yaml(args.config)
    os.makedirs(args.outdir, exist_ok=True)

    import glob
    input_files = sorted(glob.glob(os.path.join(args.input_dir, "*.tsv")))
    if not input_files:
        print(f"[cluster] No .tsv files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    for tsv_path in input_files:
        base = os.path.basename(tsv_path)
        sample = base.replace(".kept.tsv", "").replace(".parsed.tsv", "").replace(".tsv", "")
        out_csv = os.path.join(args.outdir, f"{sample}.corrected_barcodes.csv")
        rep_counts = cluster_sample(tsv_path, sample, cfg, out_csv)
        print(f"[cluster] {sample}: {len(rep_counts)} corrected barcodes, {sum(rep_counts.values())} reads", file=sys.stderr)


def cmd_map(args):
    barcode_sample_counts = load_all_counts(args.counts_dir, min_count=args.min_count)
    write_overlap_outputs(barcode_sample_counts, args.outdir)
    n_samples = len({s for c in barcode_sample_counts.values() for s in c})
    print(f"[map] {len(barcode_sample_counts)} distinct barcodes across {n_samples} samples -> {args.outdir}", file=sys.stderr)


def cmd_run(args):
    cfg_path = args.config
    extract_dir = os.path.join(args.outdir, "1_extracted")
    filter_dir = os.path.join(args.outdir, "2_filtered")
    cluster_dir = os.path.join(args.outdir, "3_clustered")
    map_dir = os.path.join(args.outdir, "4_mapped")

    extract_args = argparse.Namespace(
        config=cfg_path, samples=args.samples, fastq_dir=args.fastq_dir,
        outdir=extract_dir, mode=args.mode,
    )
    cmd_extract(extract_args)

    cluster_input_dir = extract_dir
    if args.index_map:
        filter_args = argparse.Namespace(
            config=cfg_path, parsed_dir=extract_dir, index_map=args.index_map, outdir=filter_dir,
        )
        cmd_filter(filter_args)
        cluster_input_dir = os.path.join(filter_dir, "kept")
    else:
        print("[run] No --index-map given, skipping the internal-index filter step "
              "(clustering directly on extracted barcodes).", file=sys.stderr)

    cluster_args = argparse.Namespace(config=cfg_path, input_dir=cluster_input_dir, outdir=cluster_dir)
    cmd_cluster(cluster_args)

    if args.do_map:
        map_args = argparse.Namespace(counts_dir=cluster_dir, outdir=map_dir, min_count=args.map_min_count)
        cmd_map(map_args)
    else:
        print("[run] --do-map not set, skipping cross-sample mapping step.", file=sys.stderr)


def build_parser():
    p = argparse.ArgumentParser(prog="barsqan", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    common_cfg = dict(default=None, help="Optional YAML config overriding defaults in config.py")

    p_extract = sub.add_parser("extract", help="Extract index/UMI/barcode from FASTQ files for each sample")
    p_extract.add_argument("--samples", required=True, help="Text file, one sample name per line")
    p_extract.add_argument("--fastq-dir", required=True, help="Directory containing FASTQ(.gz) files")
    p_extract.add_argument("--outdir", required=True)
    p_extract.add_argument("--mode", choices=["auto", "se", "pe"], default="auto",
                            help="Force single-end or paired-end; default auto-detects per sample")
    p_extract.add_argument("--config", **common_cfg)
    p_extract.set_defaults(func=cmd_extract)

    p_filter = sub.add_parser("filter", help="Filter parsed barcodes by expected internal index + motif")
    p_filter.add_argument("--parsed-dir", required=True, help="Directory of {sample}.parsed.tsv from `extract`")
    p_filter.add_argument("--index-map", required=True, help="CSV with columns: sample,expected_index")
    p_filter.add_argument("--outdir", required=True)
    p_filter.add_argument("--config", **common_cfg)
    p_filter.set_defaults(func=cmd_filter)

    p_cluster = sub.add_parser("cluster", help="Collapse UMIs and cluster barcode variants")
    p_cluster.add_argument("--input-dir", required=True,
                            help="Directory of .tsv files to cluster (e.g. kept/ from `filter`, or extract outdir)")
    p_cluster.add_argument("--outdir", required=True)
    p_cluster.add_argument("--config", **common_cfg)
    p_cluster.set_defaults(func=cmd_cluster)

    p_map = sub.add_parser("map", help="Map barcode occurrences across sample count tables")
    p_map.add_argument("--counts-dir", required=True, help="Directory of {sample}.corrected_barcodes.csv from `cluster`")
    p_map.add_argument("--outdir", required=True)
    p_map.add_argument("--min-count", type=int, default=2,
                        help="Minimum per-sample count for a barcode to be included (default 2)")
    p_map.set_defaults(func=cmd_map)

    p_run = sub.add_parser("run", help="Run the full pipeline: extract -> [filter] -> cluster -> [map]")
    p_run.add_argument("--samples", required=True)
    p_run.add_argument("--fastq-dir", required=True)
    p_run.add_argument("--outdir", required=True)
    p_run.add_argument("--mode", choices=["auto", "se", "pe"], default="auto")
    p_run.add_argument("--index-map", default=None, help="If given, run the filter step; else skip it")
    p_run.add_argument("--do-map", action="store_true", help="Also run the cross-sample mapping step")
    p_run.add_argument("--map-min-count", type=int, default=2)
    p_run.add_argument("--config", **common_cfg)
    p_run.set_defaults(func=cmd_run)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
