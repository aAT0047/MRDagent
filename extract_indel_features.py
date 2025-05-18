#!/usr/bin/env python3
import os
import argparse
import subprocess
import csv
import concurrent.futures
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description="Extract indel features from BAM files")
    parser.add_argument('-i', '--input_dir',    required=True,
                        help="Directory containing input .bam files")
    parser.add_argument('-o', '--output_dir',   required=True,
                        help="Directory to write per-sample indel feature CSVs and merged summary")
    parser.add_argument('-w', '--workers',      type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    parser.add_argument('-t', '--threads',      type=int, default=4,
                        help="Threads to pass to samtools when converting BAM to SAM (default: 4)")
    return parser.parse_args()

def convert_bam_to_sam(bam_path, sam_path, threads):
    """Use samtools to convert BAM to SAM."""
    cmd = ["samtools", "view", "-h", "-o", sam_path, "-@", str(threads), bam_path]
    subprocess.run(cmd, check=True)

def parse_cigar(cigar):
    """Split CIGAR string into list of (length, op) tuples."""
    ops = []
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ''
    return ops

def extract_indel_features(sam_path):
    """
    Read SAM, filter reads (MAPQ>=30, CIGAR!='*'), then count:
      - total_reads
      - insertion_events, deletion_events
      - mean insertion size, mean deletion size
    Returns a dict of features.
    """
    read_ids = set()
    ins_sizes = []
    del_sizes = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            mapq = int(fields[4])
            cigar = fields[5]
            if cigar == '*' or mapq < 30:
                continue
            read_ids.add(fields[0])
            for length, op in parse_cigar(cigar):
                if op == 'I':
                    ins_sizes.append(length)
                elif op == 'D':
                    del_sizes.append(length)
    total_reads = len(read_ids)
    total_ins = len(ins_sizes)
    total_del = len(del_sizes)
    mean_ins = sum(ins_sizes)/total_ins if total_ins else 0.0
    mean_del = sum(del_sizes)/total_del if total_del else 0.0
    mean_indels_per_read = (total_ins + total_del)/total_reads if total_reads else 0.0

    return {
        "TotalReads": total_reads,
        "InsertionCount": total_ins,
        "DeletionCount": total_del,
        "MeanInsertionSize": round(mean_ins, 2),
        "MeanDeletionSize": round(mean_del, 2),
        "MeanIndelsPerRead": round(mean_indels_per_read, 3)
    }

def process_bam(bam_path, output_dir, threads):
    """Convert BAM to SAM, extract indel features, write CSV, and clean up."""
    sample = os.path.basename(bam_path).replace('.bam','')
    sam_path = os.path.join(output_dir, f"{sample}.sam")
    csv_path = os.path.join(output_dir, f"{sample}_indel_features.csv")
    # ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    convert_bam_to_sam(bam_path, sam_path, threads)
    feats = extract_indel_features(sam_path)
    os.remove(sam_path)

    # write per-sample CSV
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sorted(feats.keys()))
        writer.writeheader()
        writer.writerow(feats)
    return csv_path

def merge_csvs(csv_paths, merged_path):
    """Merge all per-sample CSVs into a single table."""
    fieldnames = None
    with open(merged_path, 'w', newline='') as out:
        writer = None
        for p in csv_paths:
            with open(p) as f:
                reader = csv.DictReader(f)
                if writer is None:
                    fieldnames = ["SampleID"] + reader.fieldnames
                    writer = csv.DictWriter(out, fieldnames=fieldnames)
                    writer.writeheader()
                for row in reader:
                    row_out = {"SampleID": os.path.basename(p).split('_indel_features.csv')[0]}
                    row_out.update(row)
                    writer.writerow(row_out)

def main():
    args = parse_args()
    bam_files = [os.path.join(args.input_dir, f)
                 for f in os.listdir(args.input_dir) if f.endswith('.bam')]
    per_sample_dir = os.path.join(args.output_dir, "per_sample")
    os.makedirs(per_sample_dir, exist_ok=True)

    # parallel processing
    csv_paths = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as exe:
        futures = {exe.submit(process_bam, bam, per_sample_dir, args.threads): bam for bam in bam_files}
        for fut in tqdm(concurrent.futures.as_completed(futures),
                        total=len(futures), desc="Extracting indel features"):
            csv_paths.append(fut.result())

    # merge
    merged_csv = os.path.join(args.output_dir, "all_samples_indel_features.csv")
    merge_csvs(csv_paths, merged_csv)
    print(f"Done! Merged features at: {merged_csv}")

if __name__ == "__main__":
    main()
