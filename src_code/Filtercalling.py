import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from scipy.stats import norm
import subprocess
import os
import logging
from Metrics import evamain
from global_vars import ref_fasta
# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def FilterMutectCalls_objective(output_path, base_vcf, truth_vcf, **params):
    return FilterMutectCalls_objective_internal(output_path, base_vcf, truth_vcf, **params)

def FilterMutectCalls_objective_internal(output_path, base_vcf, truth_vcf, **params):
    distance_on_haplotype = int(params['distance_on_haplotype'])
    f_score_beta = params['f_score_beta']
    false_discovery_rate = params['false_discovery_rate']
    initial_threshold = params['initial_threshold']
    log_artifact_prior = params['log_artifact_prior']
    log_indel_prior = params['log_indel_prior']
    log_snv_prior = params['log_snv_prior']
    max_events_in_region = int(params['max_events_in_region'])
    min_slippage_length = int(params['min_slippage_length'])
    pcr_slippage_rate = params['pcr_slippage_rate']
    base_tar_gz = base_vcf.replace("unfiltered", "f1r2")
    tar_gz_path = os.path.join(output_path, base_tar_gz + ".tar.gz")
    artifact_priors = os.path.join(output_path, base_tar_gz + "_artifact-priors.tar.gz")
    base_vcf_path = os.path.join(output_path, f'{base_vcf}.vcf.gz')
    
    output_filename = os.path.join(
        output_path,
        f"{base_vcf}_filtered_{distance_on_haplotype}_{f_score_beta}_"
        f"{false_discovery_rate}_{initial_threshold}_"
        f"{log_artifact_prior}_{log_indel_prior}_"
        f"{log_snv_prior}_{max_events_in_region}_"
        f"{min_slippage_length}_{pcr_slippage_rate}_tmp.vcf.gz"
    )

    java_path = "/data/home/std_12/java/jdk-22/bin/java"  # Java 路径
    gatk_jar_path = "/data/home/std_12/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"  # GATK jar 路径

    command = [
        java_path, "-jar", gatk_jar_path, "FilterMutectCalls",
        '-R', ref_fasta,
        '-V', base_vcf_path,
        '-O', output_filename,
        '--distance-on-haplotype', str(distance_on_haplotype),
        '--f-score-beta', str(f_score_beta),
        '--false-discovery-rate', str(false_discovery_rate),
        '--initial-threshold', str(initial_threshold),
        '--log-artifact-prior', str(log_artifact_prior),
        '--log-indel-prior', str(log_indel_prior),
        '--log-snv-prior', str(log_snv_prior),
        '--max-events-in-region', str(max_events_in_region),
        '--min-slippage-length', str(min_slippage_length),
        '--pcr-slippage-rate', str(pcr_slippage_rate),
        "--orientation-bias-artifact-priors", artifact_priors
    ]
    # 基于基础 VCF 文件名创建日志文件路径
    log_file = os.path.join(output_path, f"{base_vcf}_filtermutectcalls.log")
    
    # 确保日志文件存在并追加写入
    with open(log_file, 'a') as log:
        # 记录运行的命令
        log.write(f'Running command: {" ".join(command)}\n')
        
        # 执行命令并将输出追加到日志文件
        # logging.info(f'Running command: {" ".join(command)}')
        subprocess.run(command, stdout=log, stderr=log, check=True)
        # subprocess.run(command,  check=True)
   
    # subprocess.run(command, check=True)

    tp, fp, fn, f1, precision, recall, TPC, FPC, FNC = evamain(output_filename, truth_vcf)
    return f1

def FilterMutectCalls_objectivelast(output_path, base_vcf, truth_vcf, **params):
    distance_on_haplotype = int(params['distance_on_haplotype'])
    f_score_beta = params['f_score_beta']
    false_discovery_rate = params['false_discovery_rate']
    initial_threshold = params['initial_threshold']
    log_artifact_prior = params['log_artifact_prior']
    log_indel_prior = params['log_indel_prior']
    log_snv_prior = params['log_snv_prior']
    max_events_in_region = int(params['max_events_in_region'])
    min_slippage_length = int(params['min_slippage_length'])
    pcr_slippage_rate = params['pcr_slippage_rate']

    split_part = os.path.basename(truth_vcf).split('.')[0]
    output_filename = os.path.join(output_path, f"{split_part}_filtered.vcf.gz")
    base_vcf_path = os.path.join(output_path, f'{base_vcf}.vcf.gz')
    base_tar_gz = base_vcf.replace("unfiltered", "f1r2")
    tar_gz_path = os.path.join(output_path, base_tar_gz + ".tar.gz")
    artifact_priors = os.path.join(output_path, base_tar_gz + "_artifact-priors.tar.gz")

    java_path = "/data/home/std_12/java/jdk-22/bin/java"  # Java 路径
    gatk_jar_path = "/data/home/std_12/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"  # GATK jar 路径

    command = [
        java_path, "-jar", gatk_jar_path, "FilterMutectCalls",
        '-R', ref_fasta,
        '-V', base_vcf_path,
        '-O', output_filename,
        '--distance-on-haplotype', str(distance_on_haplotype),
        '--f-score-beta', str(f_score_beta),
        '--false-discovery-rate', str(false_discovery_rate),
        '--initial-threshold', str(initial_threshold),
        '--log-artifact-prior', str(log_artifact_prior),
        '--log-indel-prior', str(log_indel_prior),
        '--log-snv-prior', str(log_snv_prior),
        '--max-events-in-region', str(max_events_in_region),
        '--min-slippage-length', str(min_slippage_length),
        '--pcr-slippage-rate', str(pcr_slippage_rate),
        "--orientation-bias-artifact-priors", artifact_priors
    ]

        # 基于基础 VCF 文件名创建日志文件路径
    log_file = os.path.join(output_path, f"{base_vcf}_filtermutectcalls.log")
    
    # 确保日志文件存在并追加写入
    with open(log_file, 'a') as log:
        # 记录运行的命令
        log.write(f'Running command: {" ".join(command)}\n')
        
        # 执行命令并将输出追加到日志文件
        # logging.info(f'Running command: {" ".join(command)}')
        subprocess.run(command, stdout=log, stderr=log, check=True)
        # subprocess.run(command,  check=True)
   

    tp, fp, fn, f1, precision, recall, TPC, FPC, FNC = evamain(output_filename, truth_vcf)
    return tp, fp, fn, f1, precision, recall, TPC, FPC, FNC

class BayesianOptimization:
    def __init__(self, f, pbounds, output_path, base_vcf, truth_vcf, random_state=None):
        self.f = f
        self.pbounds = pbounds
        self.output_path = output_path
        self.base_vcf = base_vcf
        self.truth_vcf = truth_vcf
        self.random_state = np.random.RandomState(random_state)
        self.keys = list(pbounds.keys())
        self.bounds = np.array([pbounds[key] for key in self.keys])
        self.X = []
        self.Y = []

    def init_points(self, n_init):
        for _ in range(n_init):
            x = [self.random_state.uniform(low, high) for low, high in self.bounds]
            self.X.append(x)
            params = dict(zip(self.keys, x))
            self.Y.append(self.f(self.output_path, self.base_vcf, self.truth_vcf, **params))
        # 添加特定的 X 值
        x_manual = [100, 1.0, 0.05, 0.1, -2.302585092994046, -16.11809565095832, -13.815510557964275, 2, 8, 0.1]
        self.X.append(x_manual)
        params_manual = dict(zip(self.keys, x_manual))
        self.Y.append(self.f(self.output_path, self.base_vcf, self.truth_vcf, **params_manual))

    def maximize(self, n_iter):
        kernel = Matern(nu=2.5, length_scale_bounds=(1e-2, 1e7))
        gp = GaussianProcessRegressor(kernel=kernel, random_state=self.random_state)
        gp.fit(np.array(self.X), np.array(self.Y))

        for _ in range(n_iter):
            x_next = self.propose_location(gp)
            params = dict(zip(self.keys, x_next))
            y_next = self.f(self.output_path, self.base_vcf, self.truth_vcf, **params)
            self.X.append(x_next)
            self.Y.append(y_next)
            gp.fit(np.array(self.X), np.array(self.Y))

    def propose_location(self, gp):
        def min_obj(X):
            return -self.expected_improvement(X, gp)

        bounds = np.array([self.pbounds[key] for key in self.keys])
        x0 = np.array([self.random_state.uniform(low, high) for low, high in bounds])
        res = minimize(min_obj, x0=x0, bounds=bounds, method='L-BFGS-B')
        return res.x

    def expected_improvement(self, X, gp, xi=0.01):
        X = X.reshape(-1, len(self.bounds))
        mu, sigma = gp.predict(X, return_std=True)
        mu_sample_opt = np.max(self.Y)

        with np.errstate(divide='warn'):
            imp = mu - mu_sample_opt - xi
            Z = imp / sigma
            ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
            ei[sigma == 0.0] = 0.0

        return ei

pbounds = {
    'distance_on_haplotype':  (10, 200),     # 原(50, 200) -> 下调下限
    'f_score_beta':           (0.5, 2.0),    # 原(0.8, 1.5) -> 扩大上下限
    'false_discovery_rate':   (0.001, 0.2),  # 原(0.01, 0.1) -> 更宽容
    'initial_threshold':      (0.005, 0.2),  # 原(0.05, 0.15) -> 放宽上下限
    'log_artifact_prior':     (-5.0, -1.5),  # 原(-3.0, -2.0) -> 覆盖更极端噪音
    'log_indel_prior':        (-20.0, -10.0),# 原(-18.0, -14.0) -> 适度扩大
    'log_snv_prior':          (-20.0, -10.0),# 原(-15.0, -12.0) -> 同理
    'max_events_in_region':   (2, 10),       # 原(2, 5) -> 允许更多事件
    'min_slippage_length':    (4, 12),       # 原(6, 10) -> 稍微扩大范围
    'pcr_slippage_rate':      (0.005, 0.3)   # 原(0.05, 0.2) -> 扩大上下限
}
# [100, 1.0,0.05,  0.1,-2.302585092994046,-16.11809565095832, -13.815510557964275, 2, 8, 0.1]
def FilterMutectmain(output_path, base_vcf, truth_vcf):
    optimizer = BayesianOptimization(f=FilterMutectCalls_objective, pbounds=pbounds, output_path=output_path, base_vcf=base_vcf, truth_vcf=truth_vcf, random_state=1)
    optimizer.init_points(n_init=200)
    optimizer.maximize(n_iter=20)

    best_idx = np.argmax(optimizer.Y)
    best_params = dict(zip(optimizer.keys, optimizer.X[best_idx]))

    tp, fp, fn, f1, precision, recall, TPC, FPC, FNC = FilterMutectCalls_objectivelast(output_path, base_vcf, truth_vcf, **best_params)
    best_score = f1
    # 创建最佳参数的字典，使用 pbounds 的键作为列名
    best_params_data = {param: best_params.get(param, np.nan) for param in pbounds.keys()}
    best_params_data['f1'] = best_score 
    # 提取文件名，不带扩展名
    file_name = os.path.splitext(os.path.splitext(os.path.basename(truth_vcf))[0])[0]
    output_file = os.path.join(output_path, f'{file_name}.csv')

    # 读取 CSV 文件的最后一行
    existing_df = pd.read_csv(output_file)
    
    # 获取最后一行的数据
    last_row = existing_df.iloc[-1].to_dict()
    
    # 合并最后一行的数据和新的最佳参数数据
    combined_data = {**last_row, **best_params_data}
    
    # 创建 DataFrame 并追加写入 CSV 文件
    combined_df = pd.DataFrame([combined_data])
    # 2. 删除最后一行
    existing_df = existing_df.drop(existing_df.index[-1])

    # 3. 保存更新后的 DataFrame 回到 CSV 文件中
    existing_df.to_csv(output_file, index=False)
    combined_df.to_csv(output_file, header=False,mode='a', index=False)
    
    print(f"数据已追加写入 {output_file}")

    return best_params, best_score, f1, precision, recall, TPC, FPC, FNC



# # ref_fasta = '/data/home/std_12/GRCH38ICGC/GRCh38_hla_decoy_ebv.fa'
# input_bam = '/data/home/std_12/ICGCCRAM/split_bam/sample_sv_27.bam'
# output_path = '/data/home/std_12/ICGCCRAM/runbam/'
# base_tar_gz = 'sample_sv_27_f1r2'
# base_vcf = 'sample_sv_27_unfiltered'
# truth_vcf_path = "/data/home/std_12/ICGCCRAM/split_vcf/sample_sv_27.vcf.gz"
#     # 生成 BAM 索引文件
# # bam_index_command = [
# #     "samtools", "index", input_bam]
# # subprocess.run(bam_index_command, check=True)
# from condiction import f1,f2
# constr_func=[f1,f2]
# # state= [0.3,0.3]
# # cleanup_files(output_path,input_bam)
# best = FilterMutectmain(output_path, base_vcf, truth_vcf_path)
# print(best)
