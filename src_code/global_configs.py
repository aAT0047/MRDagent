from threading import Event
import multiprocessing
import threading
import os
import re
import subprocess
import numpy as np
import pandas as pd
import random
from multiprocessing import Barrier,Value
from par import *
import csv

from condiction import *

synchronization_event = threading.Event()  # 用于同步任务
from  objectfun import obj_func


# 定义一个全局计数器来错开时间
counter = Value('i', 0)
# 创建全局事件对象
single_core_event = Event()
multi_core_event = Event()
checkpoint_index=False
I=1#约束采集起始
N=2 #约束数量
dim=9 #X维度
dim_T =18
n=150#初始采集样本数
K =20#ADMM迭代次数
state =[0.3,0.3]
constr_func = [f1, f2]
m=np.full(N, I*40)#每个约束的初始样本点数量
# 初始化贝叶斯优化的迭代次数
alpha = np.random.randint(1,2, K)
alpha[0] = 20# 第一次迭代赋予更多预算
beta = np.random.randint(1,2, N)
beta[0] =20 # 对于每个约束，第一次迭代赋予更多预算
rho = 0.1  # ADMM算法中的罚参数
eps = 0.01 # 收敛判定的精度
M = 1e10  # 约束罚函数中的一个系数

folder_path = '/data/home/std_12/ICGCCRAM/split_bam'
truthvcf_path = '/data/home/std_12/ICGCCRAM/split_vcf'
# output_path = '/data/home/std_12/ICGCCRAM/runbam/'
ref_fasta = '/data/home/std_12/GRCH38ICGC/GRCh38_hla_decoy_ebv.fa'
# bam_files = [f for f in os.listdir(folder_path) if f.endswith('.bam')]
# 分组参数

# # 定义全局 Barrier 对象，总任务数是
# num_tasks = len(bam_files)
# barrier = Barrier(num_tasks)
# 命令行参数解析
