
import argparse
import os
import pysam

def parse_vcf(vcf_path):
    # 初始化结果列表
    results = []
    with pysam.VariantFile(vcf_path, 'r') as vcf:
        for record in vcf:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # type: ignore # 只考虑第一个替代等位基因
            
            # 计算终止位置
            end_position = pos + len(ref) - 1 # type: ignore

            # 判断变异类型
            if len(ref) == 1 and len(alt) == 1: # type: ignore
                variant_type = '<SNV>'
            elif len(ref) > len(alt): # type: ignore
                variant_type = '<DEL>'
            else:
                variant_type = '<INS>'
            
            # 添加到结果中
            results.append((chrom, str(pos), str(end_position), variant_type))

    return results

def extract_variant_info(vcf_path):
    variant_data = []
    with pysam.VariantFile(vcf_path, 'r') as vcf:
        for record in vcf:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # type: ignore # 只考虑第一个替代等位基因
            
            # 计算终止位置
            end_position = record.stop  # END位置
            if record.info.get('SVTYPE') == 'SNV':
                variant_type = '<SNV>'
            else:
                variant_type = alt
                
            variant_data.append((chrom, str(pos), str(end_position), variant_type))
    return variant_data

def sort_vcf(input_path, output_path):
    variants = []
    with pysam.VariantFile(input_path, 'r') as vcf:
        header = vcf.header.copy()
        for record in vcf:
            variants.append(record)

    sorted_variants = sorted(variants, key=lambda x: (x.chrom, x.pos))

    with pysam.VariantFile(output_path, 'w', header=header) as out_vcf:
        for variant in sorted_variants:
            out_vcf.write(variant)

def match_score(var1, var2):
    """
    返回两个变异的匹配打分, 区间 [0.0, 1.0].
    如果完全不匹配则返回 0.0；如果非常接近则返回接近 1.0.
    具体评分逻辑可根据需求自由调整，这里只是一个示例.
    """
    _, start1, end1, type1 = var1
    _, start2, end2, type2 = var2

    start1 = int(start1)
    end1 = int(end1)
    start2 = int(start2)
    end2 = int(end2)

    # 根据类型决定允许位置误差
    if type1 == "<SNV>" or type2 == "<SNV>":
        allowed_difference = 0
    else:
        allowed_difference = 1
    # 如果类型完全不同，这里直接打 0 分（也可以酌情扣分）
    if type1 != type2:
        return 0.0

    # 计算两端位置的差值
    diff_start = abs(start1 - start2)
    diff_end = abs(end1 - end2)

    # 这里使用一个简单的线性衰减方案：
    # 当 diff == 0 时，得分最高(=1.0)；当 diff == allowed_difference 时，得分=0；中间线性插值.
    # 也可以换成更复杂的函数，比如指数衰减、阈值分段、长度占比等等。
    def position_score(diff, allowed):
        if diff >= allowed:
            return 0.0
        else:
            # 简单线性：
            return 1.0 - (diff / allowed)

    # 计算起始位点、结束位点两者的得分，可以加权平均
    score_start = position_score(diff_start, allowed_difference)
    score_end   = position_score(diff_end,   allowed_difference)

    # 这里简单取平均，也可以按自己需求加权
    final_score = 0.5 * (score_start + score_end)

    return final_score


def compute_f1(standard_vcf, called_vcf):
    """
    现在的思路是：
    - partial_TP 累加匹配分数（范围 0.0~N, N<=len(called_vcf)）
    - FP 仍然是整数，但也可以改成类似 partial_FP 的模式
    - FN 用没有得到配对(>0 分)的标准变异来计算
    """
    partial_TP = 0.0  # 用于累加打分的真阳性
    FP = 0            # 仍按原先方式记录“完全无法匹配”的次数
    matched_indices = set()  # 记录哪些 standard_vcf 的索引被匹配(部分匹配)过

    for call_var in called_vcf:
        best_score = 0.0
        best_idx = None

        # 找到在 standard_vcf 中跟 call_var 打分最高的变异
        for idx, std_var in enumerate(standard_vcf):
            score = match_score(std_var, call_var)
            if score > best_score:
                best_score = score
                best_idx = idx

        # 如果最高分 > 0，则累加到 partial_TP，并认为匹配了 best_idx
        if best_score > 0:
            partial_TP += best_score
            matched_indices.add(best_idx)
        else:
            # 否则就是彻底不匹配, FP + 1
            FP += 1

    # 没被匹配到(或者说没拿到>0分)的标准变异，都算 FN
    FN = len(standard_vcf) - len(matched_indices)

    # --- 下面是指标计算 ---
    # partial_TP 现在是一个浮点数了
    # 你可以把它当做新的 “有效匹配总量”，而 FP 依旧是个整数
    # Precision = TP / (TP + FP)，但是我们现在 TP = partial_TP
    # 也可以考虑使用 partial_TP / (partial_TP + FP)，看是否合适你的业务场景

    # 避免分母为0
    if (partial_TP + FP) == 0:
        precision = 0.0
    else:
        precision = partial_TP / (partial_TP + FP)

    if (partial_TP + FN) == 0:
        recall = 0.0
    else:
        recall = partial_TP / (partial_TP + FN)

    if (precision + recall) == 0:
        f1 = 0.0
    else:
        f1 = 2.0 * precision * recall / (precision + recall)

    # 下面几个指标也可以考虑改成基于 partial_TP 的版本
    # 这里示例保留原先“个数”概念下的 TPC, FPC, FNC
    # 你可以根据需要改写
    TPC = partial_TP / len(standard_vcf) if len(standard_vcf) > 0 else 0
    FPC = FP / len(called_vcf) if len(called_vcf) > 0 else 0
    FNC = FN / len(standard_vcf) if len(standard_vcf) > 0 else 0

    return partial_TP, FP, FN, f1, precision, recall, TPC, FPC, FNC


def evamain(input_vcf_1, input_vcf_2):

    
    sort_vcf(input_vcf_1, input_vcf_1)
    sort_vcf(input_vcf_2, input_vcf_2)
    
    # variant_info_2 = extract_variant_info(input_vcf_2)
    variant_info_2 = parse_vcf(input_vcf_2)
    variant_info_1 = parse_vcf(input_vcf_1)
    
    TP, FP, FN, f1_score, precision, recall, TPC, FPC, FNC = compute_f1(variant_info_2, variant_info_1)

    return TP, FP, FN, f1_score, precision, recall, TPC, FPC, FNC

