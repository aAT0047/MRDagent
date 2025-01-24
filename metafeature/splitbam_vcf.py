import os
import subprocess
import csv
from tqdm import tqdm

# 文件路径定义
bam_folder = "/data/home/std_12/ICGCCRAM/ctDNAMER"
vcf_folder = "/data/home/std_12/ICGCCRAM/VCF_filterbed"
bed_folder = "/data/home/std_12/ICGCCRAM/bed"
bam_output_folder = "/data/home/std_12/ICGCCRAM/split_bam"
vcf_output_folder = "/data/home/std_12/ICGCCRAM/split_vcf"
bed_output_folder = "/data/home/std_12/ICGCCRAM/split_bed"
mapping_csv = "/data/home/std_12/ICGCCRAM/mapping.csv"

# 创建输出目录
os.makedirs(bam_output_folder, exist_ok=True)
os.makedirs(vcf_output_folder, exist_ok=True)
os.makedirs(bed_output_folder, exist_ok=True)

# 初始化映射关系表
mapping_data = []

# 获取总任务数
total_tasks = sum(
    (len(open(os.path.join(bed_folder, bed_file)).readlines()) // 10)
    for bed_file in os.listdir(bed_folder)
    if bed_file.endswith(".bed")
)

# 全局计数器，用于生成全局唯一的文件名
file_counter = 1

# 使用 tqdm 进度条
with tqdm(total=total_tasks, desc="Processing", unit="task") as pbar:
    # 遍历 BED 文件
    for bed_file in os.listdir(bed_folder):
        if not bed_file.endswith(".bed"):
            continue

        sample_id = os.path.splitext(bed_file)[0]
        bed_file_path = os.path.join(bed_folder, bed_file)

        # 找到对应的 BAM 和 VCF 文件
        bam_file = os.path.join(bam_folder, f"{sample_id}.bam")
        bai_file = f"{bam_file}.bai"

        # 定义使用的核数
        threads = 4

        # 检查是否已索引
        if not os.path.exists(bai_file):
            print(f"Indexing BAM file: {bam_file}...")
            subprocess.run(["samtools", "index", "-@", str(threads), bam_file], check=True)

        vcf_file = os.path.join(vcf_folder, f"{sample_id}.vcf.gz")

        if not os.path.exists(bam_file) or not os.path.exists(vcf_file):
            print(f"BAM 或 VCF 文件缺失：{sample_id}，跳过...")
            continue

        # 按染色体划分并分割 BED 文件
        with open(bed_file_path, "r") as bed:
            bed_regions = {}
            for line in bed:
                chrom = line.split("\t")[0]
                if chrom not in bed_regions:
                    bed_regions[chrom] = []
                bed_regions[chrom].append(line)

            for chrom, lines in bed_regions.items():
                for i in range(0, len(lines), 10):
                    region_index = i // 10 + 1
                    region_file = os.path.join(bed_output_folder, f"sample_sv_{file_counter}.bed")
                    with open(region_file, "w") as region_bed:
                        region_bed.writelines(lines[i:i+10])

                    # 定义输出文件
                    bam_output_file = os.path.join(bam_output_folder, f"sample_sv_{file_counter}.bam")
                    vcf_output_file = os.path.join(vcf_output_folder, f"sample_sv_{file_counter}.vcf.gz")

                    # 使用 samtools 分割 BAM 文件并优化头文件
                    samtools_cmd = [
                        "samtools", "view", "-h", "-b", "-L", region_file, "-o", bam_output_file, bam_file,
                        "--threads", "4"
                    ]
                    subprocess.run(samtools_cmd, check=True)

                    # 为分割的 BAM 文件生成索引
                    subprocess.run(["samtools", "index", "-@", "4", bam_output_file], check=True)

                    # 使用 bcftools 分割 VCF 文件
                    bcftools_cmd = [
                        "bcftools", "view", "-R", region_file, "-o", vcf_output_file, "-Oz", vcf_file,
                        "--threads", "4"
                    ]
                    subprocess.run(bcftools_cmd, check=True)

                    # 为分割的 VCF 文件生成索引
                    subprocess.run(["bcftools", "index", "-t", vcf_output_file, "--threads", "4"], check=True)

                    # 记录分割的染色体区域信息
                    bed_region = f"{chrom}:{lines[i].split()[1]}-{lines[min(i+10, len(lines)-1)].split()[2]}"

                    # 保存映射关系
                    mapping_data.append({
                        "SampleID": sample_id,
                        "RegionIndex": region_index,
                        "NewFileName": f"sample_sv_{file_counter}",
                        "BEDRegion": bed_region,
                        "RegionFile": region_file
                    })

                    # 更新全局计数器
                    file_counter += 1

                    # 更新进度条
                    pbar.update(1)

# 保存映射关系到 CSV 文件
with open(mapping_csv, "w", newline="") as csvfile:
    fieldnames = ["SampleID", "RegionIndex", "NewFileName", "BEDRegion", "RegionFile"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(mapping_data)

print(f"所有样本分割完成，映射关系已保存到 {mapping_csv}")
