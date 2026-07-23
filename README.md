# barsqan

A single, installable package that replaces the four standalone scripts
(`extract_BC.py`, `filter_extedBC_by_intidx_BCmotif.py`,
`cluster_n_count_bcs_collapseUMI.py`, `map_barcode_overlap.py`) with one
consistent pipeline, driven by **sample names** instead of hardcoded file
paths, and by **one config object** instead of constants duplicated (and
occasionally drifting out of sync, e.g. the barcode motif) across files.

## Pipeline

```
samples.txt + FASTQ dir
        │
        ▼
 [1] extract   -- find each sample's FASTQ file(s), auto-detect SE/PE,
   │              extract Index1/UMI1/Barcode(+Index2/UMI2 for PE),
   │              validate read structure + (for PE) R1/R2 downstream
   │              overlap. Tallies unmatched ("off-target") sequences.
   ▼
 {sample}.parsed.tsv, {sample}.stats.tsv, {sample}.unmatched.tsv
        │
        ▼
 [2] filter    -- (optional) tolerant re-check: does the internal index
   │              match the *expected* index for this sample? Does the
   │              barcode still satisfy the motif with some mismatch
   │              tolerance? Flags index cross-talk / bleed.
   ▼
 kept/{sample}.kept.tsv, filtered/{sample}.filtered.tsv, filter_summary.csv
        │
        ▼
 [3] cluster   -- collapse near-identical UMIs (sequencing error), then
   │              cluster barcode sequencing-error variants into their
   │              most-abundant "true" barcode (Levenshtein distance).
   ▼
 {sample}.corrected_barcodes.csv
        │
        ▼
 [4] map       -- (optional) build a barcode x sample matrix + summary of
                  which barcodes are shared across which sample names.
   ▼
 barcode_x_sample_counts.csv, barcode_x_sample_freq.csv, barcode_overlap_summary.csv
```

Steps 1 and 3 are always run; steps 2 and 4 are optional (skip step 2 if
you don't have/need an expected-index map; skip step 4 if you only care
about per-sample barcode counts).

## Install

### Conda (recommended)

The fastest path — creates a conda environment named `barsqan`, installs all
runtime dependencies from conda-forge, and installs the `barsqan` CLI into it
in one command:

```bash
cd barsqan
conda env create -f environment.yml
conda activate barsqan
pip install --no-deps -e .
barsqan --help
```

To update the environment later (after editing `environment.yml` or pulling
new code):

```bash
conda env update -f environment.yml --prune
```

For a development environment that also includes `pytest` (and installs
`barsqan` in editable mode):

```bash
conda env create -f environment-dev.yml
conda activate barsqan-dev
pip install --no-deps -e .
pytest tests/
```

#### Building a real conda package

If you'd rather distribute `barsqan` as a conda package (so it can be
`conda install`ed like any other), a conda-build recipe is included:

```bash
conda install -n base conda-build       # if you don't already have it
conda build conda-recipe/ -c conda-forge

# then install the locally built package into a fresh env:
conda create -n barsqan -c conda-forge python">=3.9"
conda activate barsqan
conda install --use-local -c conda-forge barsqan
```

### pip

```bash
cd barsqan
pip install .          # or: pip install -e .  for an editable/dev install
```

Either route installs the `barsqan` command (and the `barsqan` Python
package, if you'd rather import the modules directly).

## Quick start

```bash
# 1) One sample name per line. Names are matched against FASTQ filenames
#    at token boundaries, so "Col-1" will NOT accidentally match "Col-10".
cat > samples.txt <<EOF
PR_pilot_Col-10_S31_L001
PR_pilot_Col-11_S32_L001
EOF

# 2) Run the whole pipeline in one command.
barsqan run \
    --samples samples.txt \
    --fastq-dir /path/to/fastqs \
    --outdir results/ \
    --index-map index_map.csv \
    --do-map \
    --config my_config.yaml
```

`index_map.csv` (only needed if you use `--index-map` / the `filter` step):

```csv
sample,expected_index
PR_pilot_Col-10_S31_L001,GCCAATTTAGGC
PR_pilot_Col-11_S32_L001,CTAAGGCATTAG
```

For single-end samples the expected index is just the 6bp `Index1` (or
whatever `index1_len` you configure); for paired-end it's `Index1+Index2`
concatenated.

Results land in:

```
results/1_extracted/   {sample}.parsed.tsv, .stats.tsv, .unmatched.tsv
results/2_filtered/    kept/{sample}.kept.tsv, filtered/{sample}.filtered.tsv, filter_summary.csv
results/3_clustered/   {sample}.corrected_barcodes.csv
results/4_mapped/      barcode_x_sample_counts.csv, barcode_x_sample_freq.csv, barcode_overlap_summary.csv
```

## Single-end vs paired-end

`--mode auto` (the default) inspects, for each sample name, how many FASTQ
files matched it:

- **two** files, one containing `_R1_`/`_R1.`/`_1.fastq` and the other
  `_R2_`/`_R2.`/`_2.fastq` → treated as **paired-end**.
- **one** file → treated as **single-end**.

You can force `--mode se` or `--mode pe` if you want extraction to fail
loudly instead of silently guessing (e.g. you expect every sample to be
paired-end and want to be told if an R2 is missing).

In single-end mode, only the R1 structure is required
(`Index1+UMI1+Primer1+Barcode`); there is no R2 downstream-overlap check
(nothing to compare against), so the pass/fail stats only report
`quality_fail` / `regex_fail`.

## Running steps individually

```bash
barsqan extract  --samples samples.txt --fastq-dir fastqs/ --outdir out/extracted --config cfg.yaml
barsqan filter   --parsed-dir out/extracted --index-map index_map.csv --outdir out/filtered --config cfg.yaml
barsqan cluster  --input-dir out/filtered/kept --outdir out/clustered --config cfg.yaml
barsqan map      --counts-dir out/clustered --outdir out/mapped
```

(If you skip `filter`, point `cluster --input-dir` at `out/extracted`
instead of `out/filtered/kept`.)

### Filtering low-abundance barcodes before mapping

Before pooling each sample's barcodes into the cross-sample overlap matrix,
you can drop barcodes whose within-sample abundance is too low *relative to
the rest of that sample* with `--min-rel-abundance`. This is a **value
threshold**: within each sample, drop every barcode whose within-sample
frequency is below `fraction × (mean|median barcode frequency)`. It removes
the long low-abundance tail (likely sequencing error / contamination) per
sample, so a barcode that is truly abundant in one sample but only a
low-abundance straggler in another is not spuriously reported as "shared":

```bash
# drop barcodes below 5% of the sample's mean barcode frequency
barsqan map --counts-dir out/clustered --outdir out/mapped \
    --min-rel-abundance 0.05 --relative-to mean --plot

# same options are available in the full pipeline
barsqan run ... --do-map --map-min-rel-abundance 0.05 --map-relative-to mean --map-plot
```

How it works and how it stacks with `--min-count`:

- `--min-count N` (default 2) is applied **first**, per barcode: any barcode
  with fewer than N reads in a sample is dropped outright.
- `--min-rel-abundance F` is then applied **per sample**: the reference
  statistic (mean or median of the sample's per-barcode within-sample
  frequencies) is computed, and every barcode whose frequency is **strictly
  below `F × reference`** is dropped. `F = 0` (default) disables it.
- `--relative-to {mean,median}` (default `median`) selects the reference
  statistic. **This choice matters a lot on skewed data:**
  - `median` is robust, but if the low-abundance tail makes up more than
    half the barcodes the median itself sits *inside* the tail, so a small
    fraction of it may drop nothing. Use a larger `F` (e.g. `0.5`) with
    median, or switch to mean.
  - `mean` is pulled up by the real high-count barcodes, so `0.05 × mean`
    typically lands in the gap between the tail and the real population —
    good for cleanly removing a large low-abundance tail.
- This is a **value threshold, not a fixed fraction**: it drops *however
  many* barcodes fall below the cutoff. That can be 5% or 90% of the
  barcodes in the file — exactly the point, if most barcodes really are
  low-abundance noise.
- A sample with only one barcode is never filtered.
- The `map` step logs the threshold and how many barcodes each sample lost:
  `[map] Col-10_S31: dropped 400/700 (57.1%) barcodes below 0.05xmean (mean freq=1.429e-03, cutoff freq=7.143e-05)`.

#### Distribution plots (`--plot`)

With `--plot` (or `--map-plot` in the full pipeline), barsqan writes one
histogram per sample to `<outdir>/plots/{sample}.count_distribution.png`.
Each plot shows the distribution of `log10(within-sample frequency)` across
that sample's barcodes, with kept barcodes (blue) and barcodes removed by
the relative-abundance filter (red) drawn as a stacked histogram. A dashed
line marks the cutoff (`F × reference`) and a dotted line marks the
reference statistic (mean/median) itself. This makes it easy to see whether
the chosen threshold cleanly separates the low-abundance tail from the real
population, or whether you should adjust `--min-rel-abundance` /
`--relative-to`.

Plotting requires `matplotlib`. It is included in the conda environment
files; if it is not installed, the `map` step still produces all CSVs and
just prints a note that plots were skipped.

## Configuring your amplicon design

Every primer sequence, index/UMI length, barcode motif, quality threshold,
overlap window, and clustering tolerance lives in one place:
[`barsqan/config.py`](barsqan/config.py). Copy
[`config.example.yaml`](config.example.yaml), edit only the keys you want
to change, and pass it with `--config` to any subcommand. Unlisted keys
keep their default.

Key ones you'll likely need to change for a new design:

| Key | Meaning |
|---|---|
| `inner_p1_primer` / `inner_p2_primer` | The primer sequences flanking the barcode / downstream region |
| `barcode_motif` | e.g. `"NNNNNacNNNNNgaNNNNNtcNNNNN"` — `N` = wildcard, lowercase = fixed base. Used to build **both** the strict extraction regex and the tolerant filter re-check, so they can't drift apart. |
| `index1_len`, `umi1_len`, `index2_len`, `umi2_len` | Fixed-length index/UMI segments |
| `qmin` | Minimum Phred quality required across index+UMI+primer(+barcode) |
| `overlap_k`, `overlap_max_mismatch` | Paired-end R1/R2 downstream-overlap validation |
| `cluster_max_dist`, `cluster_merge_ratio`, `cluster_protect_count` | Barcode error-collapsing tolerances |

## Diagnosing low pass rates

If `passed` is low relative to `total_reads` in a sample's `.stats.tsv`:

1. Check `quality_fail` vs `regex_fail` vs `overlap_fail` — they're counted
   independently and printed in that order, so you know which stage to
   investigate first.
2. Look at `{sample}.unmatched.tsv` — the most frequent R1 prefixes that
   failed the structural regex. Recurring sequences here usually point to
   an off-target amplicon, a shifted primer/index length, or a barcode
   motif that doesn't match your actual design (see the barcode-motif
   discussion in the original design conversation for an example of this
   exact failure mode).
3. For `overlap_fail`, temporarily raise `overlap_max_mismatch` or lower
   `overlap_k` to see whether reads are close-but-not-quite passing
   (likely a real design/offset issue) vs. wildly different (likely
   contamination or a wrong primer).

## Example / smoke test

`example_data/make_test_data.py` generates a small synthetic paired-end +
single-end FASTQ set (two PE samples sharing two barcodes, one SE sample)
that exercises every pipeline stage, including reads engineered to fail
each filter. Run it, then run the pipeline on the output to see the
expected shape of every output file:

```bash
cd example_data
python3 make_test_data.py
cd ..
barsqan run \
    --samples example_data/samples.txt \
    --fastq-dir example_data \
    --outdir example_data/results \
    --index-map example_data/index_map.csv \
    --do-map \
    --config config.example.yaml
```

## Troubleshooting

### Clustering seems to hang / is extremely slow

The `cluster` step corrects sequencing errors by computing edit distances
between barcodes. It uses `rapidfuzz`, which normally runs a fast compiled
C backend. On very new Python interpreters (e.g. **Python 3.14**, before
rapidfuzz publishes wheels for it) rapidfuzz silently falls back to a
**pure-Python** backend that is ~100x slower. Combined with a large number
of distinct barcodes, that can make clustering look frozen.

barsqan handles this in two ways:

- It uses a **pigeonhole block index** so each low-count barcode variant is
  only compared against the handful of candidate parents that share a
  sequence block with it — not against every other barcode. This keeps
  clustering near-linear (tens of thousands of barcodes cluster in ~1-2s
  even on the pure-Python backend).
- If the compiled rapidfuzz backend is missing, it falls back to a
  **self-contained bounded (banded) Levenshtein** with an early-exit
  cutoff, so performance never depends on a compiled extension being
  present.

The `cluster` step prints which backend is active:

```
[cluster] edit-distance backend: rapidfuzz-cpp        # fast compiled path
[cluster] edit-distance backend: pure-python-bounded  # fallback (still fast here)
```

To guarantee the fast compiled path, use a Python version rapidfuzz ships
wheels for. The provided `environment.yml` pins `python>=3.9,<3.13` for
exactly this reason. If you already built an env on Python 3.14, recreate
it:

```bash
conda deactivate
conda env remove -n barsqan
conda env create -f environment.yml   # gets a rapidfuzz-supported Python
conda activate barsqan
```

## Requirements

- Python >= 3.9 (use < 3.13 for rapidfuzz's fast compiled backend)
- `regex` (fuzzy/edit-distance regex matching)
- `rapidfuzz` (fast Levenshtein distance for barcode clustering)
- `pyyaml` (config file loading)

All three are available on both conda-forge and PyPI. Via conda they are
installed for you by `environment.yml`; via pip:

```bash
pip install -r requirements.txt
```
