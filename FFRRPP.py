import csv
import math
import re
import os
from collections import Counter
from datetime import datetime

# ================= 配置区 =================
# 输入的比对序列文件路径 (FASTA 格式)
INPUT_FASTA = "alignment.fasta"
# 输出的 DOE 清单文件名
OUTPUT_CSV = f"V02_探针优先_全基因组靶区DOE清单_{datetime.now().strftime('%Y-%m-%d')}.csv"

# 寻优参数设置
F_LEN, P_LEN, R_LEN = 20, 22, 20
MIN_GAP, MAX_GAP = 2, 25
LOCUS_WINDOW = 50  # 独立靶区物理距离阈值 (bp)
# ==========================================

def read_fasta(file_path):
    """读取并解析 FASTA 文件"""
    sequences = []
    current_seq = []
    if not os.path.exists(file_path):
        print(f"❌ 错误: 找不到文件 '{file_path}'，请确认文件路径。")
        return []
        
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append("".join(current_seq).upper())
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq).upper())
    return sequences

def calc_tm(seq):
    """计算经验 Tm (适用于长片段)"""
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    return round(64.9 + 41 * (g + c - 16.4) / len(seq), 1)

def calc_gc(seq):
    """计算 GC 含量 (%)"""
    g = seq.count('G')
    c = seq.count('C')
    return round(((g + c) / len(seq)) * 100, 1)

def is_valid_oligo(seq, is_probe=False):
    """工业级纯碱基引物综合质控"""
    if 'N' in seq or re.search(r'[^ATGC]', seq): return False
    if re.search(r'([ATGC])\1{3,}', seq): return False # 过滤 4 个以上连续重复碱基
    
    gc = calc_gc(seq)
    if gc < 40 or gc > 65: return False

    if not is_probe:
        # 引物 3' 端 GC 夹约束
        end5 = seq[-5:]
        end_gc = end5.count('G') + end5.count('C')
        if end_gc < 1 or end_gc > 3: return False
    return True

def reverse_complement(seq):
    """反向互补"""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(mapping.get(b, b) for b in reversed(seq))

def get_top_variants(start_idx, length, sequences, max_variants=2):
    """提取某窗口下真实的变异序列组合 (执行 90% 有效性硬约束)"""
    counts = Counter()
    total_valid = 0
    total_seq = len(sequences)

    for seq in sequences:
        sub = seq[start_idx:start_idx+length]
        if '-' in sub or 'N' in sub:
            continue
        counts[sub] += 1
        total_valid += 1

    # 熔断机制：有效碱基低于 90% 直接废弃该区域
    if total_valid / total_seq < 0.90:
        return []

    if not counts: return []
    
    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / total_valid

    # 混合补偿逻辑：第一大序列占比不足 97%，且第二大序列占比超过 4%
    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        second_cov = sorted_counts[1][1] / total_valid
        if second_cov > 0.04:
            variants.append(sorted_counts[1][0])
            
    return variants

def calculate_mix_mismatch(variants, start_idx, sequences):
    """联合错配统计引擎，返回详细的矩阵数据"""
    total_seq = len(sequences)
    m0 = m1 = m2 = m3p = 0
    seq_len = len(variants[0])

    for seq in sequences:
        sub = seq[start_idx:start_idx+seq_len]
        if '-' in sub or 'N' in sub:
            m3p += 1
            continue
            
        best_mismatch = seq_len
        for target_seq in variants:
            mismatch = sum(1 for a, b in zip(target_seq, sub) if a != b)
            if mismatch < best_mismatch:
                best_mismatch = mismatch
                
        if best_mismatch == 0: m0 += 1
        elif best_mismatch == 1: m1 += 1
        elif best_mismatch == 2: m2 += 1
        else: m3p += 1

    if total_seq == 0: return {}

    return {
        'p0': round((m0 / total_seq) * 100, 1),
        'p1': round((m1 / total_seq) * 100, 1),
        'p2': round((m2 / total_seq) * 100, 1),
        'p3': round((m3p / total_seq) * 100, 1),
        'm0': m0, 'm1': m1, 'm2': m2, 'm3p': m3p, 'total': total_seq
    }

def run_pipeline():
    print(f"🚀 [1/4] 正在加载比对序列库: {INPUT_FASTA}")
    sequences = read_fasta(INPUT_FASTA)
    if len(sequences) < 2:
        print("❌ 错误: 序列不足或文件格式不正确。")
        return

    seq_len = len(sequences[0])
    print(f"✅ 成功读取 {len(sequences)} 条序列，对齐长度为 {seq_len} bp。")
    print(f"🔍 [2/4] 启动滑动窗口热力学与错配矩阵引擎... (该过程可能需要几分钟)")

    all_candidates = []
    
    # 进度监控
    for i in range(seq_len - 150):
        if i % 1000 == 0 and i > 0:
            print(f"   -> 已扫描至 {i} bp...")
            
        f_variants = get_top_variants(i, F_LEN, sequences, 2)
        if not f_variants or not all(is_valid_oligo(v, False) for v in f_variants): continue

        for gap1 in range(MIN_GAP, MAX_GAP + 1):
            p_start = i + F_LEN + gap1
            p_variants = get_top_variants(p_start, P_LEN, sequences, 1)
            if not p_variants or not all(is_valid_oligo(v, True) for v in p_variants): continue

            for gap2 in range(MIN_GAP, MAX_GAP + 1):
                r_start = p_start + P_LEN + gap2
                r_variants_raw = get_top_variants(r_start, R_LEN, sequences, 2)
                if not r_variants_raw or not all(is_valid_oligo(v, False) for v in r_variants_raw): continue

                r_variants = [reverse_complement(v) for v in r_variants_raw]
                
                min_p_tm = min(calc_tm(v) for v in p_variants)
                max_f_tm = max(calc_tm(v) for v in f_variants)
                max_r_tm = max(calc_tm(v) for v in r_variants)
                
                # 热力学初筛：探针 Tm 必须高于最高引物 Tm
                if min_p_tm > max_f_tm and min_p_tm > max_r_tm: 
                    f_stats = calculate_mix_mismatch(f_variants, i, sequences)
                    p_stats = calculate_mix_mismatch(p_variants, p_start, sequences)
                    r_stats = calculate_mix_mismatch(r_variants_raw, r_start, sequences)

                    # 探针绝对优先权打分系统
                    base_score = f_stats['p0'] + (p_stats['p0'] * 3) + r_stats['p0']
                    probe_bonus = 50 if p_stats['p0'] >= 99.0 else 0
                    probe_penalty = (98.0 - p_stats['p0']) * 10 if p_stats['p0'] < 98.0 else 0
                    
                    mix_f = -8 if len(f_variants) > 1 else 0
                    mix_r = -8 if len(r_variants) > 1 else 0
                    target_p_tm = max(max_f_tm, max_r_tm) + 8
                    tm_penalty = -(abs(min_p_tm - target_p_tm) * 1.5)
                    
                    total_score = base_score + probe_bonus - probe_penalty + mix_f + mix_r + tm_penalty

                    all_candidates.append({
                        'fwd': f_variants, 'rev': r_variants, 'probe': p_variants,
                        'fStats': f_stats, 'pStats': p_stats, 'rStats': r_stats,
                        'size': r_start + R_LEN - i, 'start': i, 'score': total_score,
                    })

    print(f"📊 [3/4] 扫描完毕。提取到 {len(all_candidates)} 套初始备选方案，正在执行 NMS 空间聚类...")
    
    # 空间非极大值抑制 (NMS) 聚类
    all_candidates.sort(key=lambda x: x['score'], reverse=True)
    global_loci_groups = []
    
    for cand in all_candidates:
        found_locus = False
        for locus in global_loci_groups:
            if abs(cand['start'] - locus['anchorStart']) <= LOCUS_WINDOW:
                if len(locus['variants']) < 3: # 1 主 + 2 备
                    locus['variants'].append(cand)
                found_locus = True
                break
        if not found_locus:
            global_loci_groups.append({
                'locusId': len(global_loci_groups) + 1,
                'anchorStart': cand['start'],
                'variants': [cand]
            })

    print(f"📁 [4/4] 成功划分出 {len(global_loci_groups)} 个独立靶区。正在导出 DOE 数据矩阵...")

    # 写入 CSV 文件
    with open(OUTPUT_CSV, mode='w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(["靶区归属(Locus)", "变体角色", "综合评分", "寡核苷酸类型", "序列 (5'->3')", "Tm (°C)", "GC (%)", "完美匹配(0)", "单碱基错配(1)", "双碱基错配(2)", "废弃(≥3/缺失)", "预期产物长度 (bp)", "精确起始坐标"])
        
        for locus in global_loci_groups:
            for v_idx, cand in enumerate(locus['variants']):
                locus_name = f"靶区_{locus['locusId']}"
                role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
                score = round(cand['score'], 1)
                size = cand['size']
                start = cand['start']
                
                # 遍历并写入 Forward
                for i, seq in enumerate(cand['fwd']):
                    type_str = f"Forward_{i+1}" if len(cand['fwd']) > 1 else "Forward"
                    st = cand['fStats']
                    writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq),
                                     f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])
                
                # 遍历并写入 Probe
                for i, seq in enumerate(cand['probe']):
                    type_str = f"Probe_{i+1}" if len(cand['probe']) > 1 else "Probe"
                    st = cand['pStats']
                    writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq),
                                     f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])

                # 遍历并写入 Reverse
                for i, seq in enumerate(cand['rev']):
                    type_str = f"Reverse_{i+1}" if len(cand['rev']) > 1 else "Reverse"
                    st = cand['rStats']
                    writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq),
                                     f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])

    print(f"🎉 运行完毕！原料齐套 DOE 报告已生成至: {OUTPUT_CSV}")

if __name__ == "__main__":
    run_pipeline()