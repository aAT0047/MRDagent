#!/usr/bin/env python3
import os
import subprocess
import csv
from tqdm import tqdm
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="按 BED 区域将 BAM/VCF 分割并生成映射表"
    )
    parser.add_argument(
        "-i", "--bam_folder",
        required=True,
        help="输入 BAM 文件目录，例如 /yourpath/ctDNA"
    )
    parser.add_argument(
        "-v", "--vcf_folder",
        required=True,
        help="输入 VCF(.vcf.gz) 文件目录，例如 /yourpath/VCF"
    )
    parser.add_argument(
        "-b", "--bed_folder",
        required=True,
        help="输入 BED 文件目录，例如 /yourpath/bed"
    )
    parser.add_argument(
        "-B", "--bam_output",
        required=True,
        help="输出拆分后 BAM 文件目录，例如 /yourpath/split_bam"
    )
    parser.add_argument(
        "-V", "--vcf_output",
        required=True,
        help="输出拆分后 VCF 文件目录，例如 /yourpath/split_vcf"
    )
    parser.add_argument(
        "-D", "--bed_output",
        required=True,
        help="输出拆分后 BED 区域文件目录，例如 /yourpath/split_bed"
    )
    parser.add_argument(
        "-m", "--mapping_csv",
        required=True,
        help="映射关系 CSV 输出路径，例如 /yourpath/mapping.csv"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="samtools/bcftools 使用的线程数（可选，默认 4）"
    )
    parser.add_argument(
        "-r", "--regions_per_file",
        type=int,
        default=10,
        help="每个拆分文件包含的 BED 行数（可选，默认 10）"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # 创建输出目录
    os.makedirs(args.bam_output, exist_ok=True)
    os.makedirs(args.vcf_output, exist_ok=True)
    os.makedirs(args.bed_output, exist_ok=True)

    mapping_data = []
    file_counter = 1

    # 预先计算总任务数，用于进度条
    total_tasks = 0
    for bed_file in os.listdir(args.bed_folder):
        if bed_file.endswith(".bed"):
            path = os.path.join(args.bed_folder, bed_file)
            lines = sum(1 for _ in open(path))
            total_tasks += (lines + args.regions_per_file - 1) // args.regions_per_file

    with tqdm(total=total_tasks, desc="Processing", unit="task") as pbar:
        for bed_file in os.listdir(args.bed_folder):
            if not bed_file.endswith(".bed"):
                continue

            sample_id = os.path.splitext(bed_file)[0]
            bam_file = os.path.join(args.bam_folder, f"{sample_id}.bam")
            bai_file = bam_file + ".bai"
            vcf_file = os.path.join(args.vcf_folder, f"{sample_id}.vcf.gz")
            bed_file_path = os.path.join(args.bed_folder, bed_file)

            # 检查必需文件
            if not os.path.exists(bam_file) or not os.path.exists(vcf_file):
                print(f"[WARN] 缺少 BAM/VCF：{sample_id}，已跳过")
                continue

            # 建立 BAM 索引如果缺失
            if not os.path.exists(bai_file):
                print(f"[INFO] 索引 BAM：{bam_file}")
                subprocess.run(
                    ["samtools", "index", "-@", str(args.threads), bam_file],
                    check=True
                )

            # 读取并按染色体分组
            bed_regions = {}
            with open(bed_file_path) as bedf:
                for line in bedf:
                    chrom = line.split()[0]
                    bed_regions.setdefault(chrom, []).append(line)

            for chrom, lines in bed_regions.items():
                for i in range(0, len(lines), args.regions_per_file):
                    region_idx = i // args.regions_per_file + 1

                    # 写入子 BED 文件
                    region_bed = os.path.join(
                        args.bed_output,
                        f"{sample_id}_sv_{file_counter}.bed"
                    )
                    with open(region_bed, "w") as out_bed:
                        out_bed.writelines(lines[i:i+args.regions_per_file])

                    # 定义输出 BAM/VCF 文件
                    bam_out = os.path.join(
                        args.bam_output,
                        f"{sample_id}_sv_{file_counter}.bam"
                    )
                    vcf_out = os.path.join(
                        args.vcf_output,
                        f"{sample_id}_sv_{file_counter}.vcf.gz"
                    )

                    # samtools 分割 BAM 并索引
                    subprocess.run(
                        [
                            "samtools", "view", "-h", "-b",
                            "-L", region_bed,
                            "-o", bam_out,
                            bam_file,
                            "--threads", str(args.threads)
                        ], check=True
                    )
                    subprocess.run(
                        ["samtools", "index", "-@", str(args.threads), bam_out],
                        check=True
                    )

                    # bcftools 分割 VCF 并索引
                    subprocess.run(
                        [
                            "bcftools", "view",
                            "-R", region_bed,
                            "-o", vcf_out, "-Oz",
                            vcf_file,
                            "--threads", str(args.threads)
                        ], check=True
                    )
                    subprocess.run(
                        ["bcftools", "index", "-t", vcf_out, "--threads", str(args.threads)],
                        check=True
                    )

                    # 记录映射信息
                    start = lines[i].split()[1]
                    end   = lines[min(i + args.regions_per_file - 1, len(lines) - 1)].split()[2]
                    mapping_data.append({
                        "SampleID":    sample_id,
                        "RegionIndex": region_idx,
                        "NewFileName": f"{sample_id}_sv_{file_counter}",
                        "BEDRegion":   f"{chrom}:{start}-{end}",
                        "RegionFile":  region_bed
                    })

                    file_counter += 1
                    pbar.update(1)

    # 输出映射表
    with open(args.mapping_csv, "w", newline="") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=[
            "SampleID", "RegionIndex", "NewFileName", "BEDRegion", "RegionFile"
        ])
        writer.writeheader()
        writer.writerows(mapping_data)

    print(f"[DONE] 分割完成，映射表已写入 {args.mapping_csv}")

if __name__ == "__main__":
    main()
