import os
import subprocess
import pandas as pd
from multiprocessing import Pool

# **配置路径**
BAM_DIR = "/data/home/std_12/ICGCCRAM/split_bam"
BED_DIR = "/data/home/std_12/ICGCCRAM/split_bed"
REFERENCE_GENOME = "/data/home/std_12/GRCH38ICGC/GRCh38_hla_decoy_ebv.fa"
GERMLINE_SNP_BED = "/data/home/std_12/GRCH38ICGC/GERMLINE_SNP_GRCh38_sorted_chr.bed.gz"  # 种系 SNP BED 文件路径
THREADS = 80 # **并行进程数**

# **获取 BAM 文件列表**
bam_files = [f for f in os.listdir(BAM_DIR) if f.endswith(".bam")]

def process_bam(bam_file):
    """处理单个 BAM 文件"""
    sample_name = os.path.splitext(bam_file)[0]  # 获取样本前缀
    bam_path = os.path.join(BAM_DIR, bam_file)
    tumour_bed_path = os.path.join(BED_DIR, f"{sample_name}.bed")  # 样本专属 BED 文件

    print(f"Processing {bam_file}...")

    # 初始化变量
    ref_count, alt_count, base_quality = 0, 0, 0
    trinucleotide_counts = {}

    # # **1. 调用 samtools mpileup**
    mpileup_cmd = f"samtools mpileup -f {REFERENCE_GENOME} {bam_path}"
    mpileup_process = subprocess.run(mpileup_cmd, shell=True, capture_output=True, text=True)

    # **2. 解析 mpileup 输出**
    total_read_positions = []
    for line in mpileup_process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 6:
            # 覆盖深度
            ref_count += int(fields[3])  # 第 4 列为参考等位基因计数
            alt_count += len(fields[4])  # 第 5 列为非参考等位基因计数

            # 碱基质量
            base_quality += sum(ord(char) - 33 for char in fields[5])

            # 提取三核苷酸序列
            sequence = fields[4]
            for i in range(len(sequence) - 2):
                trinucleotide = sequence[i:i+3]
                trinucleotide_counts[trinucleotide] = trinucleotide_counts.get(trinucleotide, 0) + 1

            # 读取位置
            read_positions = [ord(char) - 33 for char in fields[5]]
            total_read_positions.extend(read_positions)

    # **3. 肿瘤样本覆盖度**
    tumour_coverage = 0
    if os.path.exists(tumour_bed_path):
        tumour_cov_cmd = f"bedtools coverage -a {tumour_bed_path} -b {bam_path}"
        tumour_coverage_output = subprocess.run(tumour_cov_cmd, shell=True, capture_output=True, text=True).stdout.strip()
        tumour_coverage = sum(int(line.split("\t")[3]) for line in tumour_coverage_output.splitlines() if line)

    # **4. 比对质量中位数**
    map_quality_cmd = f"samtools view -q 1 {bam_path} | awk '{{print $5}}'"
    map_quality_output = subprocess.run(map_quality_cmd, shell=True, capture_output=True, text=True).stdout.strip()
    map_quality_values = [int(line) for line in map_quality_output.splitlines() if line.isdigit()]
    map_quality = sorted(map_quality_values)[len(map_quality_values) // 2] if map_quality_values else 0

    # **5. 中位读长位置**
    read_position = sum(total_read_positions) / len(total_read_positions) if total_read_positions else 0

    # **6. 同聚物率 (Homopolymer Rate)**
    homopolymer_rate = 0
    if os.path.exists(tumour_bed_path):
        getfasta_cmd = f"bedtools getfasta -fi {REFERENCE_GENOME} -bed {tumour_bed_path}"
        fasta_output = subprocess.run(getfasta_cmd, shell=True, capture_output=True, text=True).stdout
        sequences = fasta_output.split("\n")
        for seq in sequences:
            if seq and not seq.startswith(">"):
                runs = [len(run) ** 2 for run in seq.split('N') if run]
                homopolymer_rate += sum(runs)
        homopolymer_rate /= len(sequences) if sequences else 1

    # **7. GC Content**
    gc_content = 0
    if os.path.exists(tumour_bed_path):
        getfasta_cmd = f"bedtools getfasta -fi {REFERENCE_GENOME} -bed {tumour_bed_path}"
        fasta_output = subprocess.run(getfasta_cmd, shell=True, capture_output=True, text=True).stdout
        gc_count = sum(seq.count("G") + seq.count("C") for seq in fasta_output.split("\n") if seq and not seq.startswith(">"))
        total_count = sum(len(seq) for seq in fasta_output.split("\n") if seq and not seq.startswith(">"))
        gc_content = gc_count / total_count if total_count else 0

    # 确保种系 SNP 文件存在
    if os.path.exists(GERMLINE_SNP_BED):
        # 运行 bedtools closest
        distance_cmd = f"bedtools closest -a {tumour_bed_path} -b {GERMLINE_SNP_BED} -d"
        # print(tumour_bed_path)
        distance_output = subprocess.run(distance_cmd, shell=True, capture_output=True, text=True).stdout.strip()
        # print(distance_output)
        # 提取距离并计算平均值
        distances = [int(line.split("\t")[-1]) for line in distance_output.splitlines() if line]
        distance_to_germline_snp = sum(distances) / len(distances) if distances else -1

    # **9. 汇总三核苷酸序列统计**
    trinucleotide_counts_str = ";".join([f"{k}:{v}" for k, v in trinucleotide_counts.items()])
    trinucleotide_total = sum(trinucleotide_counts.values())

    # 返回所有特征
    return [
        bam_file, ref_count, alt_count, base_quality, tumour_coverage, map_quality,
        read_position, homopolymer_rate, gc_content, trinucleotide_counts_str, trinucleotide_total,
        distance_to_germline_snp
    ]

# **多线程处理 BAM 文件**
with Pool(THREADS) as p:
    results = p.map(process_bam, bam_files)

# **将结果转换为 DataFrame**
df = pd.DataFrame(results, columns=[
    "BAM File", "Reference Allele Count", "Non-reference Allele Count", "Sum of Base Qualities",
    "Tumour Coverage", "Mapping Quality", "Median Read Position", "Homopolymer Rate",
    "GC Content", "Trinucleotide Sequence Counts", "Trinucleotide Total Count", "Distance to Germline SNP"
])

# **保存为 CSV 文件**
csv_path = "/data/home/std_12/ICGCCRAM/bam_features_474.csv"
df.to_csv(csv_path, index=False)
print(f"所有特征数据已保存至: {csv_path}")
