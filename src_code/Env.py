from time import sleep
import numpy as np
import gym
from gym import spaces
from MutectADMMBO  import admmmain
from condiction import f1,f2
from objectfun import*
from  FilterMutect import FilterMutectmain

import pandas as pd
import time
# from multiprocessing import 
class CustomEnv(gym.Env):
    def __init__(self, ref_fasta, input_bam, output_path, base_tar_gz, base_vcf, truth_vcf_path):
        super(CustomEnv, self).__init__()
        # 定义状态空间和动作空间
        self.observation_space = spaces.Box(low=0, high=1, shape=(2,), dtype=np.float32)
        self.action_space = spaces.Box(low=-0.1, high=0.1, shape=(2,), dtype=np.float32)
        self.ref_fasta = ref_fasta
        self.input_bam = input_bam
        self.output_path = output_path
        self.base_tar_gz = base_tar_gz
        self.base_vcf = base_vcf
        self.truth_vcf_path = truth_vcf_path
        # 初始状态
        self.state = np.array([0.3, 0.3], dtype=np.float32)
        self.iteration = 0
        self.best_params = None
        self.best_score = 0
        self.fscore = 0
        self.precision = 0
        self.recall = 0
        self.TPC = 0
        self.FPC = 0
        self.FNC = 0

    def reset(self, seed=None, options=None):
        self.state = np.array([0.3, 0.3], dtype=np.float32)
        self.iteration = 0
        print('Resetting environment...')
        return self.state, {}

    def update_metrics(self):
      
        constr_func = [f1, f2]
        self.best = admmmain(constr_func, self.state, self.ref_fasta, self.input_bam, self.output_path, self.base_tar_gz, self.base_vcf, self.truth_vcf_path)
        # print(f'{self.input_bam}ADMMBO已完成')
        self.best_params, self.best_score, self.fscore, self.precision, self.recall, self.TPC, self.FPC, self.FNC = FilterMutectmain(self.output_path, self.base_vcf, self.truth_vcf_path)
        self.F1 = self.best_score
        # print(f'{self.input_bam} Filter已完成')


    def step(self, action):
        delta_fn, delta_fp = action
        self.state[0] = np.clip(self.state[0] + delta_fn, 0, 1)
        self.state[1] = np.clip(self.state[1] + delta_fp, 0, 1)

        self.update_metrics()
        reward = self.calculate_reward(delta_fn, delta_fp)
        terminated = self.F1 > 0.8
        truncated = self.iteration >=100
        # 如果你明确需要将它们转换为布尔值，可以使用 bool 而不是 np.bool_
        terminated = bool(terminated)
        truncated = bool(truncated)
        done = bool(terminated or truncated)
        self.iteration += 1
       
        return self.state, reward, terminated, truncated, {"terminated": terminated, "truncated": truncated}

    def calculate_reward(self, delta_fn, delta_fp):
        # 计算奖励
        R_F1 = self.F1
        R_delta = -(delta_fn ** 2 + delta_fp ** 2)
        alpha = 2 * min(self.FPC, self.FNC) / (self.FPC + self.FNC) if (self.FPC + self.FNC) > 0 else 0
        R_balance = alpha

        # 权重系数
        w_1, w_2, w_3 = 1.0, 0.2, -0.3
        reward = w_1 * R_F1 + w_2 * R_delta + w_3 * R_balance
        return reward
    
    def render(self, mode='human', close=False):
        pass

# 注册自定义环境
from gym.envs.registration import register

register(
    id='CustomEnv-v0',
    entry_point='__main__:CustomEnv',
)
