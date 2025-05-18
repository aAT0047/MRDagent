#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
from multiprocessing import Pool
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract SNV features from BAM files")
    parser.add_argument(
        "-i", "--bam_dir", required=True,
        help="Directory containing input .bam files"
    )
    parser.add_argument(
        "-b", "--bed_dir", required=True,
        help="Directory containing per-sample BED files"
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Reference genome FASTA"
    )
    parser.add_argument(
        "-g", "--germline_snp_bed", required=True,
        help="Germline SNP BED (gzipped)"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=80,
        help="Number of parallel processes (default: 80)"
    )
    parser.add_argument(
        "-o", "--output_csv", required=True,
        help="Path to write the final merged CSV"
    )
    return parser.parse_args()

def process_bam(args_tuple):
    """Process a single BAM and its BED, extract SNV features."""
    bam_file, BAM_DIR, BED_DIR, REF, GERMLINE_SNP_BED = args_tuple
    sample = os.path.splitext(bam_file)[0]
    bam_path = os.path.join(BAM_DIR, bam_file)
    bed_path = os.path.join(BED_DIR, f"{sample}.bed")

    # 1. mpileup for SNV counts & base qualities
    mpileup = subprocess.run(
        ["samtools", "mpileup", "-f", REF, bam_path],
        capture_output=True, text=True, check=True
    )

    ref_count = alt_count = base_quality = 0
    trinuc_counts = {}
    read_pos_list = []

    for line in mpileup.stdout.splitlines():
        cols = line.split('\t')
        if len(cols) < 6:
            continue
        depth = int(cols[3])
        seq  = cols[4]
        qual = cols[5]

        ref_count += depth
        alt_count += len(seq)
        base_quality += sum(ord(q)-33 for q in qual)
        read_pos_list.extend(ord(q)-33 for q in qual)

        # trinucleotide contexts
        for i in range(len(seq)-2):
            tri = seq[i:i+3]
            trinuc_counts[tri] = trinuc_counts.get(tri, 0) + 1

    # 2. tumour coverage from BED + BAM
    tumour_cov = 0
    if os.path.exists(bed_path):
        cov = subprocess.run(
            ["bedtools", "coverage", "-a", bed_path, "-b", bam_path],
            capture_output=True, text=True, check=True
        ).stdout
        for l in cov.splitlines():
            parts = l.split('\t')
            tumour_cov += int(parts[3])

    # 3. mapping quality median
    mapq_out = subprocess.run(
        f"samtools view -q 1 {bam_path} | awk '{{print $5}}'",
        shell=True, capture_output=True, text=True, check=True
    ).stdout
    mapqs = [int(x) for x in mapq_out.split() if x.isdigit()]
    mapq_med = sorted(mapqs)[len(mapqs)//2] if mapqs else 0

    # 4. median read position
    read_pos_med = (
        sum(read_pos_list)/len(read_pos_list)
        if read_pos_list else 0
    )

    # 5. homopolymer rate
    hypo_rate = 0
    if os.path.exists(bed_path):
        fasta = subprocess.run(
            ["bedtools", "getfasta", "-fi", REF, "-bed", bed_path],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()
        runs = []
        for seq in fasta:
            if seq.startswith(">"): continue
            for run in seq.split('N'):
                if run:
                    runs.append(len(run)**2)
        hypo_rate = sum(runs)/len(runs) if runs else 0

    # 6. GC content
    gc = 0
    if os.path.exists(bed_path):
        fasta = subprocess.run(
            ["bedtools", "getfasta", "-fi", REF, "-bed", bed_path],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()
        seqs = [s for s in fasta if not s.startswith(">")]
        total = sum(len(s) for s in seqs)
        gc = sum(s.count("G")+s.count("C") for s in seqs) / total if total else 0

    # 7. distance to nearest germline SNP
    dist = -1
    if os.path.exists(GERMLINE_SNP_BED):
        out = subprocess.run(
            ["bedtools", "closest", "-a", bed_path, "-b", GERMLINE_SNP_BED, "-d"],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()
        ds = [int(l.split('\t')[-1]) for l in out if l]
        dist = sum(ds)/len(ds) if ds else -1

    # 8. flatten trinuc counts
    tri_str = ";".join(f"{k}:{v}" for k,v in trinuc_counts.items())
    tri_total = sum(trinuc_counts.values())

    return [
        sample, ref_count, alt_count, base_quality, tumour_cov,
        mapq_med, read_pos_med, hypo_rate, gc, tri_str, tri_total, dist
    ]

def main():
    args = parse_args()
    bam_files = [f for f in os.listdir(args.bam_dir) if f.endswith(".bam")]

    pool_args = [
        (f, args.bam_dir, args.bed_dir, args.reference, args.germline_snp_bed)
        for f in bam_files
    ]

    with Pool(args.threads) as p:
        results = list(p.map(process_bam, pool_args))

    cols = [
        "Sample", "RefCount", "AltCount", "SumBaseQ",
        "TumourCov", "MapQmed", "ReadPosMed", "HypoRate",
        "GCcontent", "TriCounts", "TriTotal", "DistToGermSNP"
    ]
    df = pd.DataFrame(results, columns=cols)
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"SNV features saved to {args.output_csv}")

if __name__ == "__main__":
    main()
