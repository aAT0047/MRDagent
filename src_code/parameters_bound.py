# 定义离散参数名列表
T_parameter_names = [
    "base-quality-score-threshold",
    "callable-depth",
    "f1r2-median-mq",
    "f1r2-min-bq",
    "max-reads-per-alignment-start",
    "pcr-indel-qual",
    "pcr-snv-qual",
    "assembly-region-padding",
    "kmer-size",
    "kmer-size",
    "max-assembly-region-size",
    "max-prob-propagation-distance",
    "min-assembly-region-size",
    "max-unpruned-variants",
    "min-dangling-branch-length",
    "phred-scaled-global-read-mismapping-rate",
    "pair-hmm-gap-continuation-penalty",
    "mbq"
]


# 定义连续参数名列表
X_parameter_names = [
    "init-lod",
    "max-af",
    "emit-lod",
    "active-probability-threshold",
    "adaptive-pruning-initial-error-rate",
    "pruning-lod-threshold",
    "flow-probability-threshold",
    "expected-mismatch-rate-for-read-disqualification",
    "min-AF"
]
# 确定 FilterMutectCalls_objectivelast 所需参数的列名
filter_parameter_names = [
    "distance_on_haplotype",
    "f_score_beta",
    "false_discovery_rate",
    "initial_threshold",
    "log_artifact_prior",
    "log_indel_prior",
    "log_snv_prior",
    "max_events_in_region",
    "min_slippage_length",
    "pcr_slippage_rate"
]

bound_T = [
    # 1. base-quality-score-threshold
    (6, 25),               # 原(6,21) -> 扩大上限到25, 允许更高阈值
    
    # 2. callable-depth
    (5, 101),              # 原(5,21) -> 可增大到100, 对超深测序更适用
    
    # 3. f1r2-median-mq
    (30, 71),              # 保持原范围, 或可下探到20 (若需要): (20,71)
    
    # 4. f1r2-min-bq
    (6, 31),               # 原(10,31) -> 下限改为6, 保留更多低质碱基
    
    # 5. max-reads-per-alignment-start
    (0, 5001),             # 原(50,5001) -> 允许“0”禁用下采样
    
    # 6. pcr-indel-qual
    [10,20,30,40,50,60],   # 原[10,20,30,40,50] -> 扩展一个60
    
    # 7. pcr-snv-qual
    [10,20,30,40,50,60],   # 同上, 扩展一个60
    
    # 8. assembly-region-padding
    (50, 2001),            # 保持原(50,2001)即可
    
    # 9. kmer-size (第一个)
    (5, 25),               # 保持或微调, 默认(5,25)足以覆盖小kmer区间
    
    # 10. kmer-size (第二个)
    (25, 50),              # 原(25,50)不变, 或可缩小到(20,50)
    
    # 11. max-assembly-region-size
    (200, 5001),           # 原(200,50001) -> 建议上限到5001, 避免过度
    
    # 12. max-prob-propagation-distance
    (40, 301),             # 原(40,1001) -> 可稍缩到(40,301)
    
    # 13. min-assembly-region-size
    (30, 151),             # 保持原(30,151)即可
    
    # 14. max-unpruned-variants
    (50, 501),             # 原(50,201) -> 上限到501, 保留更多低频变异
    
    # 15. min-dangling-branch-length
    (2, 10),               # 原(2,50) -> 无需太大; 4~8 常用
    
    # 16. phred-scaled-global-read-mismapping-rate
    (30, 51),              # 原(2,51) -> 通常30~45~60区间更常见
    
    # 17. pair-hmm-gap-continuation-penalty
    (6, 15),               # 原(6,13) -> 稍放宽到15
    
    # 18. mbq
    (6, 20)                # 原(3,20) -> 建议最低6, 防止过多噪音
]


# 连续参数边界 (bound)
# 按顺序对应 X_parameter_names
bound = [
    # 1. init-lod
    (0.5, 3.5),     # 原(1.0,3.5) -> 下限放宽到0.5, 捕捉更低频变异
    
    # 2. max-af
    (0.005, 0.05),  # 原(0.005,0.02) -> 上限扩到0.05, 防止高AF就直接过滤
    
    # 3. emit-lod
    (1.5, 3.0),     # 原(2.0,3.0) -> 适度降低下限, 提高敏感度
    
    # 4. active-probability-threshold
    (0.0005, 0.01), # 原(0.001,0.01) -> 再降低到0.0005, 避免跳过稀有变异
    
    # 5. adaptive-pruning-initial-error-rate
    (0.0005, 0.005),# 保持或可下探到(0.0001,0.005), 视测序噪音水平
    
    # 6. pruning-lod-threshold
    (2, 5),         # 保持原(2,5), 如需更敏感可减到(1,5)
    
    # 7. flow-probability-threshold
    (0.002, 0.01),  # 原值即可, IonTorrent平台若需要再微调
    
    # 8. expected-mismatch-rate-for-read-disqualification
    (0.01, 0.05) ,  # 保持原(0.01,0.05), 若质量非常高可下调到0.005

    (1e-5, 1e-3)
]



