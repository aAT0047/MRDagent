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
# import matplotlib.pyplot
from numpy import arange, meshgrid, sqrt
from objectfun import*
import random
from scipy.stats import norm
# from  global_vars import inite_dx
class problem:
  
  def __init__(self,K,dim,dim_T,n,m,delta,rho,eps,M,y1,z1,y2,z2,num_constraint,alpha,beta,bounds,bound_T,model_opt,model_feas,F,X,T,C,Z,constr_func,state,ref_fasta,input_bam,output_path,base_tar_gz,base_vcf,truth_vcf_path):
    self.ADMM_iter=K
    self.dim=dim
    self.n=n
    self.m=m
    self.delta=delta
    self.rho=rho
    self.eps=eps
    self.M=M
    self.y1=y1
    self.z1=z1
    self.num_constr=num_constraint
    self.X=X
    self.T =T
    self.beta=beta
    self.B=bounds
    self.U=[]
    self.F=F
    self.model_opt=model_opt
    self.model_feas=model_feas
    self.C=C
    self.Z=Z
    self.V =[]
    self.dim_T = dim_T
    self.B_T = bound_T
    self.y2 = y2  # T相关的拉格朗日乘子
    self.z2 =z2  # T相关的辅助变量
    self.constr_func = constr_func
    self.ref_fasta=ref_fasta
    self.input_bam=input_bam
    self.output_path=output_path
    self.base_tar_gz=base_tar_gz
    self.base_vcf=base_vcf
    self.truth_vcf_path=truth_vcf_path
    self.state = state
    self.XT =[]
 

  def initialize(self, result_Z):
      self.XT = np.hstack((self.X, self.T))
      self.U = np.zeros(self.n).reshape((self.n, -1))

        # 更新 U 值，考虑 Z 和约束，包括映射罚项和黑盒条件的约束
      for i in range(self.n):
            # 计算映射的罚项
          mapping_penalty = (self.rho / 2) * np.linalg.norm(self.T[i] - np.round(self.T[i]))**2
            # 初始化约束惩罚项
          constraint_penalty = 0
          for j in range(self.num_constr):
              constraint_penalty += (self.rho / 2) * np.linalg.norm(self.XT[i] - self.z2[j] + self.y2[j] / self.rho)**2
          # 更新 U 值
          self.U[i] = result_Z[i] + (self.rho / 2) * np.linalg.norm(self.XT[i] - self.z1 + self.y1 / self.rho)**2 + constraint_penalty + mapping_penalty

      self.initialized = True  # 标记初始化完成


  def OPT_Z(self, alpha, xstar, tstar):
    self.alpha = alpha  # 迭代次数
    # self.rho = 1.0  # 罚参数
    if not self.initialized:
        raise Exception("OPT_Z called before initialization")

    for t in range(alpha):
        # 训练高斯过程模型
        self.model_opt.fit(self.XT, self.U)

        # 生成新的 Z 样本点
        # 生成新的 X 样本点
        X_samples = np.array([np.random.uniform(low, high, 200) for (low, high) in self.B])
        # 定义每个边界的阈值比例
        threshold_ratio = 0.05
        for i, (low, high) in enumerate(self.B):
            # 计算每个边界的阈值
          threshold = (high - low) * threshold_ratio
          X_samples[i, 150:200] = np.random.uniform(max(low, xstar[i] - threshold), min(high, xstar[i] + threshold), 50)
        X_samples = np.transpose(X_samples)
        # for i, (low, high) in enumerate(self.B):
        #     X_samples[i, 150:200] = np.random.uniform(max(low, xstar[i] - 0.1), min(high, xstar[i] + 0.1), 50)
        # 

        # 生成新的 T 样本点
        T_samples = np.array([np.concatenate([np.random.uniform(b[0], b[1], 150), 
                                              np.random.uniform(max(b[0], tstar[i] - 1), min(b[1], tstar[i] + 1), 50)]) 
                              if isinstance(b, tuple) else np.random.choice(b, 200) 
                              for i, b in enumerate(self.B_T)])
        T_samples = np.transpose(T_samples)

        # 合并 X 和 T 样本点生成 Z 样本点
        Z_samples = np.hstack((X_samples, T_samples))

        # 使用高斯过程模型预测新样本点的目标函数值和标准差
        mu, std = self.surrogate_T(self.model_opt, Z_samples)
        mu = np.array(mu)
        std = np.array(std)
        best = min(self.U)

        # 计算期望改进值
        with np.errstate(divide='warn'):
            imp = mu - best - 0.01
            Z = imp / std
            ei = imp * norm.cdf(Z) + std * norm.pdf(Z)
            ei[std == 0.0] = 0.0

        # 选择具有最大期望改进的样本点
        ix = np.argmax(ei)
        selected_Z = Z_samples[ix]

        self.XT = np.vstack((self.XT, selected_Z))  # 更新 Z

        # 拆分 selected_Z 为 X 和 T
        selected_X = selected_Z[:len(xstar)]
        selected_T = selected_Z[len(xstar):]
        self.T = np.vstack((self.T, selected_T))  # 更新 T
        self.X = np.vstack((self.X, selected_X))  # 更新 X

        # 计算新的目标函数值
        updated_F,_,_ = obj_funcondiction(selected_X, selected_T, self.ref_fasta, self.input_bam, self.output_path, self.base_tar_gz, self.base_vcf, self.truth_vcf_path)

        # 更新 F
        self.F = np.vstack((self.F, updated_F))

        # 更新 U
        new_U = updated_F
        constraint_penalty =0
        mapping_penalty = (self.rho / 2) * np.linalg.norm(selected_T - np.round(selected_T))**2
        for j in range(self.num_constr):
          constraint_penalty+= (self.rho / 2) *np.linalg.norm(selected_Z - self.z2[j] + self.y2[j] / self.rho)**2
        new_U = updated_F+ (self.rho / 2) * np.linalg.norm(updated_F - self.z1 + self.y1 / self.rho)**2 + constraint_penalty + mapping_penalty                             
        self.U = np.vstack((self.U, new_U))

        # 更新 z1 和 y1
        self.z1 = selected_Z - self.y1 / self.rho
        self.y1 = self.y1 + self.rho * (selected_Z - self.z1)
        self.n += 1  # 更新样本点的计数

    # 选择并返回最佳的 Z 点
    index = np.argmin(self.U)
    best_Z = self.XT[index]
    best_X = best_Z[:len(xstar)]
    best_T = best_Z[len(xstar):]

    return best_Z, best_X,best_T,updated_F  # 返回最佳 Z 点和目标函数

  def FEAS(self, ix, x_min, t_min,inite_dx):
    H=np.zeros((self.m[ix],1))
    H.reshape((self.m[ix],-1))  
    # self.y1=y1
    # self.y2=y2
    
    x_t_min = np.concatenate((x_min, t_min))
    # print(self.m[ix])
    for i in range(self.m[ix]):
      H[i]=int(self.C[ix][i]>0)+(self.rho/(2*self.M))*(np.linalg.norm( x_t_min -self.z2[ix]+self.y2 [ix]/self.rho))**2
    hplus=min(H)
    #print(H)z
    for t in range(self.beta[ix]):
      # self.model_feas[ix].fit(self.Z[ix][0:self.m[ix]],self.C[ix][0:self.m[ix]])
      self.model_feas[ix].fit(self.Z[ix][:inite_dx[ix]-1], self.C[ix][:inite_dx[ix]-1]) 
      Zsamples,Xsamples,Tsamples=self.generate_samples(x_min, t_min)
      mu,std=self.surrogate_Z(self.model_feas[ix],Zsamples)
 
      
      probs=np.zeros(len(Zsamples))#200
      

      for i, (mu_i, std_i) in enumerate(zip(mu, std)):
        theta = 1 - norm.cdf(-mu_i / std_i)
     
        current_penalty = (self.rho / (2 * self.M)) * np.linalg.norm(Zsamples[i] - np.concatenate((x_min, t_min)) + self.y2 [ix] / self.rho) ** 2
        if hplus -current_penalty <= 0:
          probs[i]=0
        
        if hplus-current_penalty>0 and hplus-current_penalty<=1:
          probs[i]= (hplus- current_penalty )*(1-theta)
        
        if hplus-current_penalty>1:
          probs[i]= (hplus- current_penalty )*(1-theta)+theta*(hplus-1- current_penalty )

      idx=np.argmax(probs)
      #print(theta)
      # Zsamples=Zsamples.reshape(200,4)
      
      self.Z[ix][inite_dx[ix]]=Zsamples[idx]

      _,FPC,FNC = obj_funcondiction(Zsamples[idx][:self.dim],Zsamples[idx][self.dim:],self.ref_fasta,self.input_bam,self.output_path,self.base_tar_gz,self.base_vcf,self.truth_vcf_path)
      h=int(self.constr_func[ix](FNC,FPC ,self.state[ix])>0)+current_penalty

      H=np.vstack((H,np.array([h])))
      hplus=min(H)
      self.C[ix][inite_dx[ix]]=self.constr_func[ix](FNC,FPC,self.state[ix] )

      # self.m[ix]=self.m[ix]+1
      self.z2[ix] = Zsamples[idx] - self.y2[ix] / self.rho
      self.y2[ix] = self.y2[ix] + self.rho * (Zsamples[idx] - self.z2[ix])
      c = [FPC,FNC]
      inite_dx[ix] = inite_dx[ix]+1
      print(self.input_bam)
      print(inite_dx)
    zmin=np.argmin(H)
      
    zmin=(self.Z[ix][zmin])
    return self.z2[ix],self.y2[ix],c



  def generate_samples(self, x_min, t_min):
    """
    生成新的X和T的联合样本。这个方法需要根据您的具体情况实现。
    """
    # 示例实现，您需要根据实际情况调整
    # Xsamples = np.random.uniform(self.B[0][0], self.B[0][1], (200, self.dim))
    xstar = x_min
    Xsamples = np.array([np.random.uniform(low, high, 200) for (low, high) in self.B])
    threshold_ratio = 0.05
    for i, (low, high) in enumerate(self.B):
      # 计算每个边界的阈值
      threshold = (high - low) * threshold_ratio
      Xsamples[i, 50:200] = np.random.uniform(max(low, xstar[i] - threshold), min(high, xstar[i] + threshold), 150)
    Xsamples = np.transpose(Xsamples)
    # Tsamples = np.random.uniform(self.B_T[0][0], self.B_T[0][1], (200, self.dim_T))
    tstar =  t_min # 初始化T的最优值

    Tsamples_all = np.array([np.concatenate([np.random.randint(b[0], b[1]+1, 50),
                              np.random.randint(max(b[0], tstar[i]-1), min(b[1], tstar[i]+1)+1, 150)])
              if isinstance(b, tuple)
              else np.random.choice(b, 200)
              for i, b in enumerate(self.B_T)])
    Tsamples = np.transpose(Tsamples_all)

    Zsamples = np.hstack((Xsamples, Tsamples))
    return Zsamples,Xsamples,Tsamples

  def calculate_probability(self, h_delta, theta):
    """
    根据h_delta和theta计算选择当前样本的概率。
    """
    if h_delta <= 0:
        return 0
    elif 0 < h_delta <= 1:
        return h_delta * (1 - theta)
    else:
        return h_delta * (1 - theta) + theta * (h_delta - 1)

  def surrogate_X(self, model, X):
        # catch any warning generated when making a prediction
        with catch_warnings():
            # ignore generated warnings
            simplefilter("ignore")
            return model.predict(X, return_std=True)

  def surrogate_T(self, model, T):
        # catch any warning generated when making a prediction
        with catch_warnings():
            # ignore generated warnings
            simplefilter("ignore")
            return model.predict(T, return_std=True)

  def surrogate_Z(self, model, Zsamples):
        # catch any warning generated when making a prediction
        with catch_warnings():
            # ignore generated warnings
            simplefilter("ignore")
            return model.predict(Zsamples, return_std=True)

  def concat(self):
        for i in range(self.num_constr):
            self.X = np.concatenate((self.X, self.Z[i][:, :self.dim]), axis=0)
            self.T = np.concatenate((self.T, self.Z[i][:, self.dim:]), axis=0)
            # for j in range(len(self.Z[i])):
            #   self.F = np.vstack((self.F, obj_func(self.Z[i][:, :self.dim][j], self.Z[i][:, self.dim:][j])))







  