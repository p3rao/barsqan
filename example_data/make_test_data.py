#!/usr/bin/env python3
"""Generate small synthetic FASTQ files to smoke-test the bctoolkit pipeline.

Creates:
  Sample_A_S1_L001_R1_001.fastq.gz / _R2_001.fastq.gz   (paired-end)
  Sample_B_S2_L001_R1_001.fastq.gz / _R2_001.fastq.gz   (paired-end, shares 2 barcodes with Sample_A)
  Sample_C_S3_L001_R1_001.fastq.gz                       (single-end)
"""
import gzip
import os
import random

random.seed(7)

INDEX1 = "GCCAAT"
INDEX2 = "TTAGGC"
INNER_P1_PRIMER = "GCCACTCGAGCACCTAGGAGG"
INNER_P2_PRIMER = "AATCAGTCCAGTTATGCTGTGAA"
DOWNSTREAM_CONST1 = (
    "GCCACGGCCGCGTCGACCCCACGCCCCTCTTTAATACGACGGGCAATTTGCACTTCAGAAAATGAAGAGTTTGCTTTAGCCATAACAAA"
)


def rev_comp(seq):
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def random_bases(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def make_barcode():
    return random_bases(5) + "AC" + random_bases(5) + "GA" + random_bases(5) + "TC" + random_bases(5)


def make_umi():
    return random_bases(8)


HIGH_QUAL = "I"  # Phred 40


def qual(seq, char=HIGH_QUAL):
    return char * len(seq)


def write_fastq_gz(path, records):
    with gzip.open(path, "wt") as f:
        for i, (seq, q) in enumerate(records):
            f.write(f"@read{i}/1\n{seq}\n+\n{q}\n")


def build_pe_pair(barcode, umi1=None, umi2=None, bad_quality=False, bad_regex=False, bad_overlap=False):
    umi1 = umi1 or make_umi()
    umi2 = umi2 or make_umi()
    r1 = INDEX1 + umi1 + INNER_P1_PRIMER + barcode + DOWNSTREAM_CONST1
    r2 = INDEX2 + umi2 + INNER_P2_PRIMER + rev_comp(DOWNSTREAM_CONST1)

    if bad_regex:
        # scramble the primer so R1 no longer matches the expected structure
        r1 = INDEX1 + umi1 + "TTTTTTTTTTTTTTTTTTTTT" + barcode + DOWNSTREAM_CONST1

    if bad_overlap:
        # replace R2 downstream with unrelated sequence -> no shared k-mer overlap with R1 downstream
        r2 = INDEX2 + umi2 + INNER_P2_PRIMER + random_bases(len(DOWNSTREAM_CONST1))

    q1 = qual(r1)
    q2 = qual(r2)
    if bad_quality:
        qc_len = 6 + 8 + len(INNER_P1_PRIMER) + 26
        q1 = "#" * qc_len + qual(r1[qc_len:])  # '#' = Phred 2, well below qmin

    return (r1, q1), (r2, q2)


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))

    shared_bc_1 = make_barcode()
    shared_bc_2 = make_barcode()
    sampleA_only_bc = make_barcode()
    sampleB_only_bc = make_barcode()

    # ---------------- Sample_A (paired-end) ----------------
    a_r1, a_r2 = [], []
    for _ in range(30):
        (r1, q1), (r2, q2) = build_pe_pair(shared_bc_1)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))
    for _ in range(15):
        (r1, q1), (r2, q2) = build_pe_pair(shared_bc_2)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))
    for _ in range(10):
        (r1, q1), (r2, q2) = build_pe_pair(sampleA_only_bc)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))
    # a few reads that should fail each filter stage, to exercise the counters
    for _ in range(5):
        (r1, q1), (r2, q2) = build_pe_pair(make_barcode(), bad_quality=True)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))
    for _ in range(5):
        (r1, q1), (r2, q2) = build_pe_pair(make_barcode(), bad_regex=True)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))
    for _ in range(5):
        (r1, q1), (r2, q2) = build_pe_pair(make_barcode(), bad_overlap=True)
        a_r1.append((r1, q1)); a_r2.append((r2, q2))

    write_fastq_gz(os.path.join(outdir, "Sample_A_S1_L001_R1_001.fastq.gz"), a_r1)
    write_fastq_gz(os.path.join(outdir, "Sample_A_S1_L001_R2_001.fastq.gz"), a_r2)

    # ---------------- Sample_B (paired-end, shares barcodes with A) ----------------
    b_r1, b_r2 = [], []
    for _ in range(8):
        (r1, q1), (r2, q2) = build_pe_pair(shared_bc_1)
        b_r1.append((r1, q1)); b_r2.append((r2, q2))
    for _ in range(20):
        (r1, q1), (r2, q2) = build_pe_pair(shared_bc_2)
        b_r1.append((r1, q1)); b_r2.append((r2, q2))
    for _ in range(12):
        (r1, q1), (r2, q2) = build_pe_pair(sampleB_only_bc)
        b_r1.append((r1, q1)); b_r2.append((r2, q2))

    write_fastq_gz(os.path.join(outdir, "Sample_B_S2_L001_R1_001.fastq.gz"), b_r1)
    write_fastq_gz(os.path.join(outdir, "Sample_B_S2_L001_R2_001.fastq.gz"), b_r2)

    # ---------------- Sample_C (single-end) ----------------
    c_r1 = []
    se_bc = make_barcode()
    for _ in range(20):
        umi1 = make_umi()
        r1 = INDEX1 + umi1 + INNER_P1_PRIMER + se_bc + DOWNSTREAM_CONST1
        c_r1.append((r1, qual(r1)))
    write_fastq_gz(os.path.join(outdir, "Sample_C_S3_L001_R1_001.fastq.gz"), c_r1)

    print("Wrote test FASTQ files to", outdir)
    print("shared_bc_1 =", shared_bc_1)
    print("shared_bc_2 =", shared_bc_2)
    print("sampleA_only_bc =", sampleA_only_bc)
    print("sampleB_only_bc =", sampleB_only_bc)
    print("se_bc (Sample_C) =", se_bc)


if __name__ == "__main__":
    main()
