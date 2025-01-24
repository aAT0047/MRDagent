from math import sin
from math import pi
import numpy as np
from numpy import arange
from numpy import vstack
from numpy import argmax
from numpy import asarray
from numpy.random import normal
from numpy.random import random
import pandas as pd
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from warnings import catch_warnings
from warnings import simplefilter
import math
from objectfun import*
import random
from problem import *
from global_vars import*
from par import bound_T,bound
import os
import time
import csv

def read_samples_and_results_from_csv(file_path):
    prob_C = np.zeros((N, np.max(m) + np.max(beta) * K))
    prob_Z = np.zeros((N, np.max(m) + K * np.max(beta), dim + dim_T))
    X_initial = np.zeros((n, dim))
    T_initial = np.zeros((n, dim_T))
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    prob_C_section = []
    result_Z_section = []
    X_initial_section = []
    T_initial_section = []
    prob_Z_section = []
    current_section = None
    skip_line = False
    
    for line in lines:
        line = line.strip()
        if line.startswith("# prob_C"):
            current_section = "prob_C"
            skip_line = True
        elif line.startswith("# prob_Z"):
            current_section = "prob_Z"
            skip_line = True
        elif line.startswith("# result_Z"):
            current_section = "result_Z"
            skip_line = True
        elif line.startswith("# X_initial"):
            current_section = "X_initial"
            skip_line = True
        elif line.startswith("# T_initial"):
            current_section = "T_initial"
            skip_line = True
        else:
            if skip_line:
                skip_line = False
                continue
            if current_section == "prob_C":
                prob_C_section.append(line)
            elif current_section == "prob_Z":
                prob_Z_section.append(line)
            elif current_section == "result_Z":
                result_Z_section.append(line)
            elif current_section == "X_initial":
                X_initial_section.append(line)
            elif current_section == "T_initial":
                T_initial_section.append(line)
    
    # 使用 pandas 读取数据，并转换为 NumPy 数组
    prob_C = pd.DataFrame([x.split(',') for x in prob_C_section], dtype=float).to_numpy().reshape((N, np.max(m) + np.max(beta) * K))
    
    result_Z = pd.DataFrame([x.split(',') for x in result_Z_section], dtype=float).to_numpy()
    X_initial = pd.DataFrame([x.split(',') for x in X_initial_section], dtype=float).to_numpy().reshape((n, dim))
    T_initial = pd.DataFrame([x.split(',') for x in T_initial_section], dtype=float).to_numpy().reshape((n, dim_T))
    prob_Z = pd.DataFrame([x.split(',') for x in prob_Z_section], dtype=float).to_numpy().reshape((N, np.max(m) + K * np.max(beta), dim + dim_T))
    return prob_C, prob_Z, result_Z, X_initial, T_initial

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')
def conditionsamples(prob, bound, bound_T):

    x_samples = np.empty((prob.num_constr, max(prob.m), prob.dim))
    t_samples = np.empty((prob.num_constr, max(prob.m), prob.dim_T))
    
    for i in range(prob.num_constr):
        for j in range(prob.m[i]):
            # 为每个约束生成 x_sample
            x_sample = np.array([np.random.uniform(low, high) for low, high in bound])
            x_samples[i][j] = x_sample

            # 为每个约束生成 t_sample
            t_sample = np.array([
                random.choice([x for x in range(b[0], b[1]+1, 10)]) if isinstance(b, tuple) and b == (30, 50) 
                else random.choice(b) if isinstance(b, list) 
                else random.randint(b[0], b[1]) 
                for b in bound_T
            ])
            t_samples[i][j] = t_sample

            # 将 x_sample 和 t_sample 分配给 Z
            prob.Z[i][j][:prob.dim] = x_sample
            prob.Z[i][j][prob.dim:] = t_sample

    return  prob.Z, x_samples, t_samples

def process_conditionsamples(prob, bound, bound_T, obj_func, constr_func, state, ref_fasta, input_bam, output_path, base_tar_gz, base_vcf, truth_vcf_path,I):
    # 生成采样点
    prob.Z, x_samples, t_samples = conditionsamples(prob, bound, bound_T)
    
    # 调用 obj_func 并行处理
    _, FPC_list, FNC_list = obj_func(
        x_samples.reshape(-1, prob.dim),
        t_samples.reshape(-1, prob.dim_T),
        ref_fasta,
        input_bam,
        output_path,
        base_tar_gz,
        base_vcf,
        truth_vcf_path,
        I
    )
    # 将 FPC 和 FNC 转换为 NumPy 数组并重塑回与 prob.C 匹配的形状
    FPC = np.array(FPC_list).reshape(prob.num_constr, max(prob.m))
    FNC = np.array(FNC_list).reshape(prob.num_constr, max(prob.m))
    # time.sleep(700)  # 在每个任务之间添加10分钟的间隔

    # 将 FPC 和 FNC 按顺序输入到 prob.C
    for i in range(prob.num_constr):
        for j in range(prob.m[i]):
            prob.C[i][j] = constr_func[i](FPC[i][j], FNC[i][j],state[i])
    
    return prob.C



def ADMMBO(dim,dim_T, n, N, K, bound,bound_T,constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path,I):
    # 序列化任务队列
    # 初始化参数
    ref_fasta=ref_fasta
    input_bam=input_bam
    output_path=output_path
    base_tar_gz=base_tar_gz
    base_vcf=base_vcf
    truth_vcf_path=truth_vcf_path
    dim=dim
    dim_T = dim_T
    n=n
    state =state
    I=I
    # m = np.full(N, I*4)  # 每个约束的初始样本点数量
   
    bound=bound
    bound_T =bound_T
  
    # 创建初始的拉格朗日乘子和对偶变量
    # 初始化 y 和 z

    y1 =np.full((N, dim + dim_T), 0.0)
    y2 =np.full((N, dim + dim_T), 0.0)
    z1 =np.full((N, dim + dim_T), 0.0)
    z2 =np.full((N, dim + dim_T), 0.0)

    
    delta=0.05
    # 初始化高斯过程回归模型
    model_opt = GaussianProcessRegressor()
    model_feas = [GaussianProcessRegressor() for i in range(N)]
    
    # 初始化目标函数值数组和设计变量数组
    F=np.zeros(n)
    F=F.reshape((n,-1))
    X = np.zeros((n,dim))
    T = np.zeros((n,dim_T))
    C=np.zeros((N, np.max(m)+np.max(beta)*K))
  
    Z = np.zeros((N, np.max(m)+K*np.max(beta), dim+dim_T))
    
    # 创建problem实例
    prob = problem(K,dim,dim_T,n,m,delta,rho,eps,M,y1,z1,y2,z2,N,alpha,beta,bound,bound_T,model_opt,model_feas,F,X,T,C,Z,constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path)


  
    csv_file_path = os.path.join(output_path, f"{os.path.basename(input_bam)}_initial_samples.csv")
    prob.C,prob.Z ,result_Z ,prob.X,prob.T= read_samples_and_results_from_csv(csv_file_path)


    # 通知所有任务 obj_func 调用已完成
    synchronization_event.set()
    synchronization_event.wait()
    synchronization_event.clear()
    # 初始采样点的目标函数值
    answers = np.zeros((prob.ADMM_iter + 1, prob.dim + prob.dim_T))  # 注意这里维度的调整以适应X和T
    
    # 初始化解和对偶残差
    r = np.full((N, 1), 1e9)  # 假设有N个约束，r初始设置为高值
    s = np.full((N, 1), 1e9)  # s同上

    mini = 1e9 + 0.1  # 用于跟踪最小目标函数值的变量
    zprev_X = np.zeros(z1.shape)  # 存储上一轮迭代中关于X的辅助变量的值
    zprev_T = np.zeros(z2.shape)  # 存储上一轮迭代中关于T的辅助变量的 # 存储上一轮迭代中关于T的辅助变量的值
    zprev_Z = np.zeros(z1.shape+z2.shape)
    xstar = np.array([2.0,0.01,3.0,0.002,0.001,2.302585092994046,0.003,0.02,0.0])  # 初始化X的最优值
    tstar =  np.array([18,10,50,20,50,40,40,100,10,25,300,50,50,100,4,45,10,10])  # 初始化T的最优值
    S = False  # 标记是否找到解

    # 修改迭代过程，修改约束函数增加Tindex
    # ADMM迭代过程
    k = 0
    answers = []
 
    synchronization_event.set()
    synchronization_event.wait()
    synchronization_event.clear()
    prob.initialize(result_Z)
    # print('初始化完成')
    # mutect2(ref_fasta, input_bam,f'{output_path}{base_tar_gz}_default.tar.gz',f'{output_path}{base_vcf}_default.vcf.gz')
    # _, _, _, _, _, _,defaultTPC,_,_= evamain(f'{output_path}{base_vcf}_default.vcf.gz', truth_vcf_path)
    
    while k < prob.ADMM_iter and not S:
        # # 1. 针对 X 执行最优化子问题
        # xstar = prob.OPT_X(alpha[k], y1, z1, xstar, tstar)
        # clear_screen()
        print(f"Current value of k: {input_bam}_{k}")
        # # 2. 针对 T 执行最优化子问题
        # tstar ,current_f= prob.OPT_T(alpha[k], y2, z2, tstar, xstar)
        best_Z, xstar, tstar, current_f = prob.OPT_Z(alpha[k],xstar, tstar)
        inite_dx= [prob.m[0],prob.m[0]]
        for i in range(prob.num_constr):
            # 3. 同时更新 X 和 T 相关的对偶变量和拉格朗日乘子
            # 注意：这里的 FEAS 函数需要能够同时处理 X 和 T
            z2[i], y2[i] ,c= prob.FEAS(i,xstar, tstar,inite_dx)
            # 4. 计算原始和对偶残差
            r[i] = max(0,c[i]-prob.state[i]) 
            s[i]  = prob.rho *np.linalg.norm(z2[i] -   zprev_Z[i])
            zprev_Z =z2[i]
            # file_path = '/home/cloudam/DQNadmm/results.txt'
            # with open( file_path, 'a') as f:
            #     f.write(f"Iteration {k}, Index {i}:\n")
            #     f.write(f"r[i]: {r[i]}\n")
            #     f.write(f"s[i]: {s[i]}\n")
            #     # f.write(f"SS[i]: {SS[i]}\n")
            #     f.write("\n")

        # 5. 检查是否满足收敛条件
        if np.all(np.linalg.norm(r) <= eps) :
          S = True  # 收敛标志设置为True
          current_f=  obj_funcsingel(xstar, tstar,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path)
          best ={'xstar': xstar, 'tstar': tstar,'current_f': current_f}
          break
        # 6. 动态调整罚参数
        if np.linalg.norm(r) >= np.linalg.norm(s) * 10:
            prob.rho *= 2
        elif np.linalg.norm(s) >= np.linalg.norm(r) * 10:
            prob.rho /= 2
        
        # 7. 记录当前迭代的解
            # 记录解和目标函数值
        result_dict = {'xstar': xstar, 'tstar': tstar, 'current_f': current_f}
        answers.append(result_dict)
        k += 1  # 迭代次数增加
        # 循环结束后查找最小的目标函数值及其对应的解
    if not S:
        # 找到所有距离 0 最近的解
        min_distance = min(abs(answer['current_f']) for answer in answers)
        closest_answers = [answer for answer in answers if abs(answer['current_f']) == min_distance]

        # 尝试排除初始化值作为解
        excluded_xstar = np.array([2.0, 0.01, 3.0, 0.002, 0.001, 2.302585092994046, 0.003, 0.02, 0.0])
        excluded_tstar = np.array([18, 10, 50, 20, 50, 40, 40, 100, 10, 25, 300, 50, 50, 100, 4, 45, 10, 10])
        
        # 筛选非默认解
        non_default_answers = [
            answer for answer in closest_answers
            if not (np.array_equal(answer['xstar'], excluded_xstar) and 
                    np.array_equal(np.round(answer['tstar']), excluded_tstar))
        ]

        # 如果有非默认解，选择其中一个；否则使用默认解
        best_solution = non_default_answers[0] if non_default_answers else closest_answers[0]

        # 构造最优解
        best = {
            'xstar': best_solution['xstar'],
            'tstar': np.round(best_solution['tstar']),
            'current_f': best_solution['current_f']
        }

        # 调用目标函数验证解
        current_f = obj_funcsingel(
            best['xstar'], 
            np.round(best['tstar']),
            ref_fasta, 
            input_bam, 
            output_path, 
            base_tar_gz, 
            base_vcf, 
            truth_vcf_path
        )
    return best

    


def admmmain(constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path):
    ref_fasta=ref_fasta
    input_bam=input_bam
    output_path=output_path
    base_tar_gz=base_tar_gz
    base_vcf=base_vcf
    truth_vcf_path=truth_vcf_path
    state = state
    best = ADMMBO(dim,dim_T,n,N,K,bound,bound_T,constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path,I)
    print (best)
    # 合并参数名称和值
    parameters = X_parameter_names + T_parameter_names
    values = np.concatenate((best['xstar'], best['tstar']))
    # 创建字典以合并所有数据
    data = {param: value for param, value in zip(parameters, values)}
    data['current_f'] = best['current_f']
    data['input_bam'] = input_bam
    
    # 创建 DataFrame
    best_df = pd.DataFrame([data])
    # 调整列顺序，将 'input_bam' 列放到第一列
    columns = ['input_bam'] + [col for col in best_df.columns if col != 'input_bam']
    best_df = best_df[columns]
    print( best_df)
    # 生成文件名
    # 提取文件名，不带扩展名
    file_name = os.path.splitext(os.path.splitext(os.path.basename(truth_vcf_path))[0])[0]
    output_file = os.path.join(output_path, f'{file_name}.csv')
    
    best_df.to_csv(output_file, mode='a', header=False, index=False)

    print(f"数据已写入 {output_file}")

    return best


# ref_fasta = '/home/cloudam/my_folder_graph/ref.fasta'
# input_bam = '/data/home/std_12/ICGCCRAM/split_bam/sample_sv_27.bam'
# output_path = '/data/home/std_12/ICGCCRAM/runbam/'
# base_tar_gz = 'sample_sv_27_f1r2'
# base_vcf = 'sample_sv_27_unfiltered'
# truth_vcf_path = "/data/home/std_12/ICGCCRAM/split_vcf/sample_sv_27.vcf"
#     # 生成 BAM 索引文件
# bam_index_command = [
#     "samtools", "index", input_bam]
# subprocess.run(bam_index_command, check=True)
# from condiction import f1,f2
# constr_func=[f1,f2]
# state= [0.3,0.3]
# # cleanup_files(output_path,input_bam)
# best = admmmain(constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path)
# print(best)
