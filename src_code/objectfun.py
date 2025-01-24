import subprocess
from math import sin
from math import pi
from numpy import arange
from numpy import vstack
from numpy import argmax
from numpy import asarray
from numpy.random import normal
from numpy.random import random
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from warnings import catch_warnings
from warnings import simplefilter
# from matplotlib import pyplot
import numpy as np
# import matplotlib.pyplot as plt
import math
import multiprocessing
from numpy import arange, meshgrid, sqrt
# from evaf1 import*
from Metrics import evamain
from par import X_parameter_names,T_parameter_names
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import time
import pysam
# from global_vars import output_path
from collections import Counter

def get_most_common_chromosome(input_bam):
    """
    通过统计 BAM 文件中 reads 的染色体分布，返回 reads 数最多的染色体。
    """
    try:
        # 使用 samtools 提取 BAM 文件的所有 reads，并统计染色体分布
        result = subprocess.run(
            ["samtools", "view", input_bam],
            capture_output=True,
            text=True,
            check=True
        )
        chromosomes = Counter()
        for line in result.stdout.splitlines():
            chrom = line.split("\t")[2]  # 第 3 列为染色体名称
            chromosomes[chrom] += 1
        # 找出 reads 数量最多的染色体
        most_common_chromosome, _ = chromosomes.most_common(1)[0]
        return most_common_chromosome
    except subprocess.CalledProcessError as e:
        print(f"Error reading BAM file: {e}")
        return None


def process_vcfs(params_list, truth_vcf_path):
    """
    Process each VCF, evaluate it, and collect TP values.
    """
    tpC_values = []
    FNC_values = []
    FPC_values = []
    
    for index, params in enumerate(params_list):
        output_vcf_gz = params['output_vcf']
        
        if not os.path.exists(output_vcf_gz):
            print(f"File {output_vcf_gz} not found. Skipping this file.")
            tpC_values.append(0.5)
            FNC_values.append(0.5)
            FPC_values.append(0.5)
            continue
        
        try:
            tp, fp, fn, f1, precision, recall, TPC, FPC, FNC = evamain(output_vcf_gz, truth_vcf_path)
            tpC_values.append(1 - TPC)
            FNC_values.append(FNC)
            FPC_values.append(FPC)
            # print(f"Processed {output_vcf_gz}: TP={tp}, FP={fp}, FN={fn}, F1={f1}, Precision={precision}, Recall={recall}")
            # print(f'成功获得第 {index + 1} 采集点')
        except Exception as e:
            print(f"Error processing {output_vcf_gz}: {e}")
            tpC_values.append(0.5)
            FNC_values.append(0.5)
            FPC_values.append(0.5)
    
    return tpC_values, FNC_values, FPC_values

def unzip_vcf(vcf_gz_path, output_dir):
    """
    Unzip a gzipped VCF file to a specified output directory.
    """
    command = ['gunzip', '-c', vcf_gz_path]
    output_path = os.path.join(output_dir, os.path.basename(vcf_gz_path)[:-3])
    with open(output_path, 'wb') as f_out:
        subprocess.run(command, stdout=f_out, check=True)
    return output_path
def mutect2(ref_fasta, input_bam, tar_gz, output_vcf, **kwargs):
    """
    Run GATK Mutect2 with specified parameters and additional optional parameters.
    
    Parameters:
    - ref_fasta: Path to the reference fasta file.
    - input_bam: Path to the input BAM file.
    - tar_gz: Path for the output tar.gz file.
    - output_vcf: Path for the output VCF file.
    - kwargs: Additional command line arguments for GATK Mutect2.
    """
    java_path = "/data/home/std_12/java/jdk-22/bin/java"  # Java 路径
    gatk_jar_path = "/data/home/std_12/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"  # GATK jar 路径
    
    # 构建命令
    command = [
        java_path, "-jar", gatk_jar_path, "Mutect2",
        "-R", ref_fasta,
        "-I", input_bam,
        "--f1r2-tar-gz", tar_gz,
        "-O", output_vcf,
        "--cloud-prefetch-buffer", "80",
        "--cloud-index-prefetch-buffer", "80",
        "--native-pair-hmm-threads", "16",  # 添加逗号
        # "--enable-dynamic-read-disqualification", "True",  # 修改为正确的参数格式
        # "--disable-tool-default-read-filters", "True"  # 修改为正确的参数格式
    ]

    # Optional: Check if additional-param is needed in each call
    # if 'additional_param' not in kwargs:
    #     command.extend(["--additional-param", "value"])
    
    # Add additional optional parameters from kwargs
    for key, value in kwargs.items():
        if isinstance(value, list):
            # If the value is a list, add each element separately
            for val in value:
                command.extend([f"--{key.replace('_', '-')}", str(val)])
        elif value is not None:
            # If the value is not a list, add it directly
            command.extend([f"--{key.replace('_', '-')}", str(value)])
    # Execute the command
        # 基于输入 BAM 文件名创建日志文件路径
    bam_basename = os.path.basename(input_bam)
    bam_filename = os.path.splitext(bam_basename)[0]
    output_path = output_vcf[:output_vcf.rfind('/') + 1]
    log_file = os.path.join(output_path,f"{bam_filename}_mutect2.log")

    
    # 打开日志文件进行追加写入
    with open(log_file, 'a') as log:
        # 执行命令并将输出追加到日志文件
        subprocess.run(command, stdout=log, stderr=log, check=True)
        # subprocess.run(command, check=True)

def run_mutect2(params):
    """
    Wrapper function to execute mutect2 with provided parameters.
    """
    mutect2(**params)

def run_task_with_event(params):
    # while not multi_core_event.is_set():
    #     if single_core_event.is_set():
    #         # print("Single-core phase active. Task is paused.")
    #         time.sleep(10)  # 暂停任务，直到 multi_core_event 被设置
    #     else:
    #         break  # multi_core_event 被设置，跳出循环
    run_mutect2(params)

def obj_func(X, T, ref_fasta, input_bam, output_path, base_tar_gz, base_vcf, truth_vcf_path,I):
    """
    Main function to set up and execute multiple Mutect2 tasks in parallel.
    """
    kwargs_list = []
    T = T.astype(int)
    # 循环遍历每个 X 和 T 的元素
    for i in range(len(X)):
        X_parameters_dict = {name: value for name, value in zip(X_parameter_names, X[i])}
        # T_parameters_dict = {name: value for name, value in zip(T_parameter_names, T[i])}
        T_parameters_dict = {}; [(T_parameters_dict[name].append(value) if isinstance(T_parameters_dict.get(name), list) else T_parameters_dict.update({name: [T_parameters_dict[name], value]}) if name in T_parameters_dict else T_parameters_dict.update({name: value})) for name, value in zip(T_parameter_names, T[i])]
        # print(T_parameters_dict)


        # 合并两个字典
        kwarg = {**X_parameters_dict, **T_parameters_dict}
        
        # 添加到 kwargs 列表
        kwargs_list.append(kwarg)

    params_list = []
    for index, kwargs in enumerate(kwargs_list, start=I):
        tar_gz = os.path.join(output_path, f'{base_tar_gz}_{index}_tmp.tar.gz')
        vcf_gz = os.path.join(output_path, f'{base_vcf}_{index}_tmp.vcf.gz')
        chromosome =  get_most_common_chromosome(input_bam)
        params = {
            'ref_fasta': ref_fasta,
            'input_bam': input_bam,
            'tar_gz': tar_gz,
            'output_vcf':  vcf_gz,
            'L': chromosome ,
            **kwargs
        }
        params_list.append(params)

    num_cores = os.cpu_count()
    with ProcessPoolExecutor(max_workers=num_cores) as executor: # type: ignore
        futures = [executor.submit(run_task_with_event, params) for params in params_list]
        results = []
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(f"Exception during parallel execution: {e}")

    eva_list = [
        {'output_vcf': f'{output_path}{base_vcf}_{i+1}_tmp.vcf.gz'}
        for i in range(len(kwargs_list))
    ]

    tpC_list, fnC_list, fpC_list = process_vcfs(eva_list, truth_vcf_path)
    # cleanup_files(output_path, input_bam)
    # time.sleep(100) 
    return tpC_list, fnC_list, fpC_list
#  obj_func是并行的再写一个单核的

def run_learn_read_orientation_model(tar_gz_path):
    """
    Run LearnReadOrientationModel on the given tar.gz file.

    Parameters:
    - tar_gz_path: Path to the input tar.gz file.

    Returns:
    - Path to the generated artifact priors file.
    """
    java_path = "/data/home/std_12/java/jdk-22/bin/java"  # Java 路径
    gatk_jar_path = "/data/home/std_12/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"  # GATK jar 路径
    # 确保输出路径以 .tar.gz 结尾
    artifact_priors = tar_gz_path.replace('.tar.gz', '_artifact-priors.tar.gz')

    command = [
        java_path, "-jar", gatk_jar_path, "LearnReadOrientationModel",
        "-I", tar_gz_path,
        "-O", artifact_priors
    ]

    print(f"[INFO] Running LearnReadOrientationModel on {tar_gz_path}")
    subprocess.run(command, check=True)

    return artifact_priors


def obj_funcsingel(X,T,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path):
   
    # 循环遍历每个 X 和 T 的元素
    T = T.astype(int)
    X_parameters_dict = {name: value for name, value in zip(X_parameter_names, X)}
    T_parameters_dict = {name: value for name, value in zip(T_parameter_names, T)}
    chromosome =  get_most_common_chromosome(input_bam)
    kwargs = {**X_parameters_dict, **T_parameters_dict}
        

    params = {
            'ref_fasta': ref_fasta,
            'input_bam': input_bam,
            'tar_gz': f'{output_path}{base_tar_gz}.tar.gz',
            'output_vcf': f'{output_path}{base_vcf}.vcf.gz',
            'L': chromosome ,
            **kwargs
    }
    

    run_mutect2(params)
    artifact_priors = run_learn_read_orientation_model(f'{output_path}{base_tar_gz}.tar.gz')
    # output_vcf = unzip_vcf(f'{output_path}{base_vcf}.vcf.gz',output_path)
    tp, fp, fn, f1, precision, recall ,TPC,FPC,FNC= evamain(f'{output_path}{base_vcf}.vcf.gz', truth_vcf_path)

    return 1-TPC

#  obj_func是并行的再写一个单核的
import hashlib

def obj_funcondiction(X, T, ref_fasta, input_bam, output_path, base_tar_gz, base_vcf, truth_vcf_path):
    # 将 T 转换为整数
    T = T.astype(int)
    
    # 创建参数字典
    X_parameters_dict = {name: value for name, value in zip(X_parameter_names, X)}
    T_parameters_dict = {name: value for name, value in zip(T_parameter_names, T)}
    
    # 合并两个字典
    kwargs = {**X_parameters_dict, **T_parameters_dict}
    
    # 生成唯一的文件标识符
    unique_id = hashlib.md5(str(kwargs).encode()).hexdigest()
    chromosome =  get_most_common_chromosome(input_bam)
    # 创建参数字典
    params = {
        'ref_fasta': ref_fasta,
        'input_bam': input_bam,
        'tar_gz': f'{output_path}{base_tar_gz}_{unique_id}_tmp.tar.gz',
        'output_vcf': f'{output_path}{base_vcf}_{unique_id}_tmp.vcf.gz',
        'L': chromosome ,
        **kwargs
    }
    
    # 运行 Mutect2
    run_mutect2(params)
    
    # 计算评价指标
    tp, fp, fn, f1, precision, recall, TPC, FPC, FNC = evamain(f'{output_path}{base_vcf}_{unique_id}_tmp.vcf.gz', truth_vcf_path)
    
    # 清理临时文件
    # cleanup_files(params_list)  # 如果有清理函数，请调用
    
    return 1-TPC,FPC, FNC




# # 示例使用
# input_bam = "/data/home/std_12/ICGCCRAM/split_bam/sample_sv_48.bam"

# ref_fasta = '/data/home/std_12/GRCH38ICGC/GRCh38_hla_decoy_ebv.fa'

# output_path = '/data/home/std_12/ICGCCRAM/runbam/'
# base_tar_gz = 'sample_sv_48_f1r2'
# base_vcf = 'sample_sv_48_unfiltered'
# truth_vcf_path = "/data/home/std_12/ICGCCRAM/split_vcf/sample_sv_48.vcf.gz"
# X= np.array([2.0, 0.01, 3.0, 0.002, 0.001, 2.302585092994046, 0.003, 0.02,0.0])
# T = np.array([18, 10, 50, 20, 50, 40, 40, 100, 10, 25, 300, 50, 50, 100, 4, 45, 10, 10])
# tt = obj_funcsingel(X,T,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path)
# print(tt )