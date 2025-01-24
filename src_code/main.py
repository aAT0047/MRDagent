import subprocess
import sys
from DQN import *
import os
import re
import concurrent.futures
import psutil
import time
import asyncio
import logging
# from global_vars import single_core_event, multi_core_event,synchronization_event
from tqdm import tqdm
from rich.console import Console
from rich.progress import Progress
from rich.logging import RichHandler

from global_vars import *
from objectfun import*
import subprocess
import random
import pandas as pd
import argparse
# 配置 rich 日志
console = Console()
logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(console=console)])

# # 命令行参数解析
parser = argparse.ArgumentParser(description="Process BAM files for a specific group.")
parser.add_argument("--bam_files_path", type=str, required=True, help="Path to the TXT file containing the BAM files list")
parser.add_argument("--output_path", type=str, required=True, help="Path to save the results for this group")
args = parser.parse_args()

# 提取参数
bam_files_path = args.bam_files_path
output_path =args.output_path

# 确保输出路径存在
os.makedirs(output_path, exist_ok=True)

# 读取 BAM 文件列表
with open(bam_files_path, 'r') as f:
    bam_files = [line.strip() for line in f if line.strip()]  # 去除空行

print(f"BAM files to process: {bam_files}")
print(f"Output directory: {output_path}")




# 确保输出路径存在
os.makedirs(output_path, exist_ok=True)

import os

def retain_csv_files(output_path):
    """
    Retain only .csv files in the specified directory. All other files will be deleted.
    """
    # 检查路径是否存在
    if not os.path.exists(output_path):
        print(f"Directory does not exist: {output_path}")
        return

    for file_name in os.listdir(output_path):
        file_path = os.path.join(output_path, file_name)
        # 检查文件是否以 .csv 结尾
        if not file_name.endswith('.csv'):
            try:
                # 删除非 .csv 文件
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Removed: {file_path}")
                else:
                    print(f"File does not exist: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")
        else:
            print(f"Retained: {file_path}")  # 保留 .csv 文件的日志

def index_bam_file(bam_file):
    # 构造 BAM 索引文件路径
    bai_file = bam_file + ".bai"
    
    # 检查 .bai 文件是否存在
    if os.path.exists(bai_file):
        print(f"Index file already exists: {bai_file}")
    else:
        print(f"Index file not found: {bai_file}. Creating index...")
        index_command = ['samtools', 'index', bam_file]
        subprocess.run(index_command, check=True)
    
def generate_initial_samples(bam_file, I):
    split_number = re.search(r'sample_sv_(\d+)\.bam', bam_file).group(1)  # type: ignore # 提取分割编号
    input_bam = os.path.join(folder_path, bam_file)
    index_bam_file(input_bam)
    base_tar_gz = f'sample_sv_{split_number}bam_f1r2'
    base_vcf = f'sample_sv_{split_number}bam_unfiltered'
    truth_vcf_path = os.path.join(truthvcf_path , f'sample_sv_{split_number}.vcf.gz')
    
    x_samples = np.empty((N, max(m), dim))
    t_samples = np.empty((N, max(m), dim_T))
    prob_Z = np.zeros((N, np.max(m) + K * np.max(beta), dim + dim_T))
    
    for i in range(N):
        for j in range(m[i]):
            x_sample = np.array([np.random.uniform(low, high) for low, high in bound])
            x_samples[i][j] = x_sample

            t_sample = np.array([
                random.choice([x for x in range(b[0], b[1] + 1, 10)]) if isinstance(b, tuple) and b == (30, 50)
                else random.choice(b) if isinstance(b, list)
                else random.randint(b[0], b[1])
                for b in bound_T
            ])
            t_samples[i][j] = t_sample

            prob_Z[i][j][:dim] = x_sample
            prob_Z[i][j][dim:] = t_sample
    # 初始化X和T的最优值
    xstar = np.array([2.0, 0.01, 3.0, 0.002, 0.001, 2.302585092994046, 0.003, 0.02,0.0])
    tstar = np.array([18, 10, 50, 20, 50, 40, 40, 100, 10, 25, 300, 50, 50, 100, 4, 45, 10, 10])
    # 随机选择一个点进行替换
    rand_i = random.randint(0, N - 1)
    rand_j = random.randint(0, m[rand_i] - 1)

    # 替换 x_samples 和 t_samples 中的随机点
    x_samples[rand_i][rand_j] = xstar
    t_samples[rand_i][rand_j] = tstar

    # 更新 prob_Z 中对应的位置
    prob_Z[rand_i][rand_j][:dim] = xstar
    prob_Z[rand_i][rand_j][dim:] = tstar

    _, FNC_list, FPC_list = obj_func(
        x_samples.reshape(-1, dim),
        t_samples.reshape(-1, dim_T),
        ref_fasta,
        input_bam,
        output_path,
        base_tar_gz,
        base_vcf,
        truth_vcf_path,
        I
    )

    prob_C = np.zeros((N, np.max(m) + np.max(beta) * K))
    
    FPC = np.array(FPC_list).reshape(N, max(m))
    FNC = np.array(FNC_list).reshape(N, max(m))
    X_initial = np.zeros((n, dim))
    T_initial = np.zeros((n, dim_T))
    
    for i in range(N):
        for j in range(m[i]):
            prob_C[i][j] = constr_func[i](FPC[i][j], FNC[i][j], state[i])
    
    for i in range(n):
        X_initial[i] = np.array([np.random.uniform(low, high) for low, high in bound])
        T_initial[i] = np.array([
            random.choice([x for x in range(b[0], b[1] + 1, 10)]) if isinstance(b, tuple) and b == (30, 50)
            else random.choice(b) if isinstance(b, list)
            else random.randint(b[0], b[1])
            for b in bound_T
        ])
    # 替换 X_initial 和 T_initial 中的随机点
    X_initial[rand_i] = xstar
    T_initial[rand_i] = tstar

    index = N * m[0] + 1

    result_Z, _, _ = obj_func(X_initial, T_initial, ref_fasta, input_bam, output_path, base_tar_gz, base_vcf, truth_vcf_path, index)
    result_Z = np.array(result_Z)  # 确保 result_Z 是 NumPy 数组

    # 将三维数组 prob_Z 展平为二维数组
    prob_Z_flat = prob_Z.reshape(-1, prob_Z.shape[-1])

    # 将结果写入CSV文件
    csv_file_path = os.path.join(output_path, f"{os.path.basename(input_bam)}_initial_samples.csv")

    dfprob_C = pd.DataFrame(prob_C)
    dfprob_Z = pd.DataFrame(prob_Z_flat)  # 展平后的三维数组
    dfresult_Z = pd.DataFrame(result_Z)
    dfX_initial = pd.DataFrame(X_initial)
    dfT_initial = pd.DataFrame(T_initial)
    headers = ['# prob_C', '# prob_Z', '# result_Z', '# X_initial', '# T_initial']

    with open(csv_file_path, 'w', newline='') as f:
        for df, header in zip([dfprob_C, dfprob_Z, dfresult_Z, dfX_initial, dfT_initial], headers):
            f.write(f'{header}\n')
            df.to_csv(f, index=False, header=True)
    
    return csv_file_path




def find_best_rows(output_path):
    all_best_rows = []

    # 获取所有 CSV 文件列表
    csv_files = [f for f in os.listdir(output_path) if f.endswith('.csv')]

    # 使用 tqdm 创建红色进度条
    for filename in tqdm(csv_files, desc="Processing files", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}', colour='red'):
        file_path = os.path.join(output_path, filename)
        
        try:
            # 读取 CSV 文件
            df = pd.read_csv(file_path)
            
            # 检查 'f1' 列是否存在并且不是空的
            if 'f1' in df.columns and not df['f1'].empty:
                # 找到 F1 最大值的那一行
                best_row = df.loc[df['f1'].idxmax()]
                all_best_rows.append(best_row)
            else:
                print(f"Skipping file {filename}: 'f1' column is missing or empty.")
        except Exception as e:
            print(f"Skipping file {filename}: {e}")

    # 合并所有最佳行到一个 DataFrame
    # print(type(all_best_rows), all_best_rows)

    all_best_rows_df = pd.DataFrame(all_best_rows)
    return all_best_rows_df

def export_best_rows_to_csv(output_path):
    all_best_rows_df = find_best_rows(output_path)
    output_file = os.path.join(output_path, 'final.csv')
    # 写入 CSV 文件
    all_best_rows_df.to_csv(output_file, index=False)
    print(f"数据已写入 {output_file}")

# 清理临时文件
def cleanup_tmp_files(output_path):
    """
    Remove all files in the specified directory that include '_tmp' in their name.
    """
    for file_name in os.listdir(output_path):
        if '_tmp' in file_name:
            file_path = os.path.join(output_path, file_name)
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Removed: {file_path}")



def process_bam_file(bam_file):
    # cleanup_tmp_files(output_path)
    split_number =re.search(r'sample_sv_(\d+)\.bam', bam_file).group(1)
    input_bam = os.path.join(folder_path, bam_file)
    # bam_index_command = ["samtools", "index", input_bam]
    # subprocess.run(bam_index_command, check=True)
    base_tar_gz = f'sample_sv_{split_number}bam_f1r2'
    base_vcf = f'sample_sv_{split_number}bam_unfiltered'
    truth_vcf_path = os.path.join(truthvcf_path , f'sample_sv_{split_number}.vcf.gz')
    env_kwargs = {
        'ref_fasta': ref_fasta,
        'input_bam': input_bam,
        'output_path': output_path,
        'base_tar_gz': base_tar_gz,
        'base_vcf': base_vcf,
        'truth_vcf_path': truth_vcf_path
    }
    # 训练模型并保存最佳状态
    best_model_state = train_dqn(1, env_kwargs)
    torch.save(best_model_state, f'{base_vcf}best_model.pth')
    logging.info(f"Trained model for {bam_file} with best state saved.")
    logging.info(f"文件 {bam_file} 已完成.")
cpu_count = psutil.cpu_count()

async def monitor_resources(futures, progress, task):
    while not all(future.done() for future in futures):
        current_cpu_usage = psutil.cpu_percent(interval=1)
        current_memory_usage = psutil.virtual_memory().percent
        logging.info(f"Current CPU usage: {current_cpu_usage}%, Memory usage: {current_memory_usage}%")

        if current_cpu_usage < 10 and current_memory_usage < 10:
            single_core_event.set()
            multi_core_event.clear()
            # logging.info("Single-core phase active. Pausing tasks.")
        else:
            single_core_event.clear()
            multi_core_event.set()
            logging.info("Multi-core phase active. Allowing tasks.")

        # 更新进度条
        done_count = sum(future.done() for future in futures)
        progress.update(task, completed=done_count)
        # await asyncio.sleep(1)  # 小休息，避免阻塞事件循环
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
def main():
    single_core_event.clear()
    multi_core_event.set()
    synchronization_event.clear()  # 重置同步事件
    logging.info("Multi-core phase active. Allowing tasks.")

   
    print('初始化开始')
    for bam_file in tqdm(bam_files, desc="Processing BAM files"):
        csv_file_path = generate_initial_samples(bam_file, I)
 
    print('初始化已完成')


    with concurrent.futures.ThreadPoolExecutor(max_workers=45) as executor:
        futures = [executor.submit(process_bam_file, bam_file) for bam_file in bam_files]
        
        with Progress(console=console) as progress:
            task = progress.add_task("Processing BAM files", total=len(futures))
            
            while not all(future.done() for future in futures):
                current_cpu_usage = psutil.cpu_percent(interval=1)
                current_memory_usage = psutil.virtual_memory().percent
                # 更新进度条并打印当前进度
                done_count = sum(future.done() for future in futures)
                progress.update(task, completed=done_count)
                print(f"Processing BAM files: {done_count}/{len(futures)} completed.")

                
                # 更新进度条
                done_count = sum(future.done() for future in futures)
                progress.update(task, completed=done_count)
            for future in futures:
                result = future.result()
headers = [
    "bam","init-lod", "max-af", "emit-lod", "active-probability-threshold", "adaptive-pruning-initial-error-rate", 
    "pruning-lod-threshold", "flow-probability-threshold", "expected-mismatch-rate-for-read-disqualification", "min-AF",
    "base-quality-score-threshold", "callable-depth", "f1r2-median-mq", "f1r2-min-bq", 
    "max-reads-per-alignment-start", "pcr-indel-qual", "pcr-snv-qual", "assembly-region-padding", 
    "kmer-size", "max-assembly-region-size", "max-prob-propagation-distance", "min-assembly-region-size", 
    "max-unpruned-variants", "min-dangling-branch-length", "phred-scaled-global-read-mismapping-rate", 
    "pair-hmm-gap-continuation-penalty", "mbq", "current_f", "distance_on_haplotype", "f_score_beta", 
    "false_discovery_rate", "initial_threshold", "log_artifact_prior", "log_indel_prior", 
    "log_snv_prior", "max_events_in_region", "min_slippage_length", "pcr_slippage_rate", "f1"
]
# 遍历每个 BAM 文件并生成相应的 CSV 文件
for bam_file in bam_files:
    # 构建相应的 CSV 文件名
    csv_file = bam_file.replace('.bam', '.csv')
    
    # CSV 文件的路径
    csv_path = os.path.join(output_path, csv_file)
    
    # 检查文件是否存在，如果不存在则写入标题
    if not os.path.isfile(csv_path):
        with open(csv_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
        print(f'Created CSV file with headers: {csv_path}')
    else:
        print(f'CSV file already exists: {csv_path}')

if __name__ == "__main__":
    main()
    print('生成参数策略文件final.csv。。。。')
    export_best_rows_to_csv(output_path)
    cleanup_tmp_files(output_path)
    retain_csv_files(output_path)
    print('一切都结束了，晚安玛卡巴卡')