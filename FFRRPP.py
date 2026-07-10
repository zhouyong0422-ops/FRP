import streamlit as st
import csv
import re
import io
from collections import Counter
from datetime import datetime

# ================= 核心参数配置 =================
F_LEN, P_LEN, R_LEN = 20, 22, 20
MIN_GAP, MAX_GAP = 2, 25
LOCUS_WINDOW = 50  # 独立靶区物理距离阈值 (bp)
# ================================================

# --- 页面基础设置 ---
st.set_page_config(page_title="V02 原料齐套 - 智能引探寻优系统", layout="wide")

def read_fasta_from_string(fasta_string):
    """从文本流中解析 FASTA 序列"""
    sequences = []
    current_seq = []
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                sequences.append("".join(current_seq).upper())
                current_seq = []
        elif line:
            current_seq.append(line)
    if current_seq:
        sequences.append("".join(current_seq).upper())
    return sequences

def calc_tm(seq):
    """经验 Tm 公式 (长片段适用)"""
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    return round(64.9 + 41 * (g + c - 16.4) / len(seq), 1)

def calc_gc(seq):
    """GC 含量计算"""
    g = seq.count('G')
    c = seq.count('C')
    return round(((g + c) / len(seq)) * 100, 1)

def is_valid_oligo(seq, is_probe=False):
    """工业级纯碱基综合质控"""
    if 'N' in seq or re.search(r'[^ATGC]', seq): return False
    if re.search(r'([ATGC])\1{3,}', seq): return False 
    
    gc = calc_gc(seq)
    if gc < 40 or gc > 65: return False

    if not is_probe:
        end5 = seq[-5:]
        end_gc = end5.count('G') + end5.count('C')
        if end_gc < 1 or end_gc > 3: return False
    return True

def reverse_complement(seq):
    """反向互补互换"""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(mapping.get(b, b) for b in reversed(seq))

def get_top_variants(start_idx, length, sequences, max_variants=2):
    """90% 覆盖率硬熔断及混合变体提取"""
    counts = Counter()
    total_valid = 0
    total_seq = len(sequences)

    for seq in sequences:
        sub = seq[start_idx:start_idx+length]
        if '-' in sub or 'N' in sub:
            continue
        counts[sub] += 1
        total_valid += 1

    if total_seq == 0 or total_valid / total_seq < 0.90:
        return []

    if not counts: return []
    
    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / total_valid

    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        second_cov = sorted_counts[1][1] / total_valid
        if second_cov > 0.04:
            variants.append(sorted_counts[1][0])
            
    return variants

def calculate_mix_mismatch(variants, start_idx, sequences):
    """多变体联合错配矩阵统计"""
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

def build_csv_string(global_loci_groups):
    """构建用于下载的 CSV 内存流"""
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["靶区归属(Locus)", "变体角色", "综合评分", "寡核苷酸类型", "序列 (5'->3')", "Tm (°C)", "GC (%)", "完美匹配(0)", "单碱基错配(1)", "双碱基错配(2)", "废弃(≥3/缺失)", "预期产物长度 (bp)", "精确起始坐标"])
    
    for locus in global_loci_groups:
        for v_idx, cand in enumerate(locus['variants']):
            locus_name = f"靶区_{locus['locusId']}"
            role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
            score = round(cand['score'], 1)
            size = cand['size']
            start = cand['start']
            
            for i, seq in enumerate(cand['fwd']):
                type_str = f"Forward_{i+1}" if len(cand['fwd']) > 1 else "Forward"
                st_data = cand['fStats']
                writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq), f"{st_data['p0']}% ({st_data['m0']}/{st_data['total']})", f"{st_data['p1']}% ({st_data['m1']}/{st_data['total']})", f"{st_data['p2']}% ({st_data['m2']}/{st_data['total']})", f"{st_data['p3']}% ({st_data['m3p']}/{st_data['total']})", size, start])
            
            for i, seq in enumerate(cand['probe']):
                type_str = f"Probe_{i+1}" if len(cand['probe']) > 1 else "Probe"
                st_data = cand['pStats']
                writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq), f"{st_data['p0']}% ({st_data['m0']}/{st_data['total']})", f"{st_data['p1']}% ({st_data['m1']}/{st_data['total']})", f"{st_data['p2']}% ({st_data['m2']}/{st_data['total']})", f"{st_data['p3']}% ({st_data['m3p']}/{st_data['total']})", size, start])

            for i, seq in enumerate(cand['rev']):
                type_str = f"Reverse_{i+1}" if len(cand['rev']) > 1 else "Reverse"
                st_data = cand['rStats']
                writer.writerow([locus_name, role, score, type_str, seq, calc_tm(seq), calc_gc(seq), f"{st_data['p0']}% ({st_data['m0']}/{st_data['total']})", f"{st_data['p1']}% ({st_data['m1']}/{st_data['total']})", f"{st_data['p2']}% ({st_data['m2']}/{st_data['total']})", f"{st_data['p3']}% ({st_data['m3p']}/{st_data['total']})", size, start])
    return output.getvalue()


# ================= Web 界面与主循环 =================
st.title("🧬 自动化引探寻优与覆盖率验证面板")
st.markdown("**(V9.0 探针绝对优先 | 90% 缺失熔断 | NMS 空间聚类 | 组合引物补偿)**")

uploaded_file = st.file_uploader("📂 导入对齐后的 FASTA 序列库", type=['fasta', 'fas', 'txt', 'aln'])

if uploaded_file is not None:
    fasta_string = uploaded_file.getvalue().decode("utf-8")
    sequences = read_fasta_from_string(fasta_string)
    
    if len(sequences) < 2:
        st.error("❌ 序列读取失败，或文件包含的序列少于2条，请检查文件格式。")
    else:
        seq_len = len(sequences[0])
        st.success(f"✅ 成功读取 {len(sequences)} 条序列，对齐总长度：{seq_len} bp。")
        
        if st.button("🚀 执行工业级全景寻优", type="primary"):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            all_candidates = []
            total_steps = seq_len - 150
            
            status_text.text("🔍 正在启动热力学与错配扫描引擎...")
            
            for i in range(total_steps):
                # UI 进度更新
                if i % 50 == 0 or i == total_steps - 1:
                    progress = int((i / total_steps) * 100)
                    progress_bar.progress(progress)
                    status_text.text(f"🔍 扫描进度：正在评估位点 {i} / {total_steps} ...")
                    
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
                        
                        if min_p_tm > max_f_tm and min_p_tm > max_r_tm: 
                            f_stats = calculate_mix_mismatch(f_variants, i, sequences)
                            p_stats = calculate_mix_mismatch(p_variants, p_start, sequences)
                            r_stats = calculate_mix_mismatch(r_variants_raw, r_start, sequences)

                            # 评分系统
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
                                # 修复了这里的拼写错误：把 -probePenalty 改成了 -probe_penalty
                                'details': {'base': base_score, 'pBonus': probe_bonus, 'pPenalty': -probe_penalty, 'mixF': mix_f, 'mixR': mix_r, 'tmPenalty': tm_penalty}
                            })

            progress_bar.empty()
            status_text.text("⚙️ 扫描结束，正在执行 NMS 空间非极大值抑制聚类...")
            
            # 空间非极大值抑制 (NMS)
            all_candidates.sort(key=lambda x: x['score'], reverse=True)
            global_loci_groups = []
            
            for cand in all_candidates:
                found_locus = False
                for locus in global_loci_groups:
                    if abs(cand['start'] - locus['anchorStart']) <= LOCUS_WINDOW:
                        if len(locus['variants']) < 3:
                            locus['variants'].append(cand)
                        found_locus = True
                        break
                if not found_locus:
                    global_loci_groups.append({
                        'locusId': len(global_loci_groups) + 1,
                        'anchorStart': cand['start'],
                        'variants': [cand]
                    })

            if not global_loci_groups:
                st.error("⚠️ 体系设计失败。在严格硬约束下，未能找到完美的扩增靶区。")
            else:
                st.success(f"🎉 寻优完成！共提取出 {len(global_loci_groups)} 个独立黄金靶区。")
                
                # 构建 CSV 下载文件 (使用 utf-8-sig 防止 Excel 乱码)
                csv_str = build_csv_string(global_loci_groups)
                st.download_button(
                    label="📥 一键导出完整 DOE 清单 (Excel CSV)",
                    data=csv_str.encode('utf-8-sig'),
                    file_name=f"V02_全基因组靶区DOE清单_{datetime.now().strftime('%Y-%m-%d')}.csv",
                    mime="text/csv",
                    type="primary"
                )

                st.markdown("---")
                st.markdown("### 🎯 黄金靶区概览 (Top 5 展示)")
                
                # 前端预览展示前 5 个 Locus
                for locus in global_loci_groups[:5]:
                    with st.expander(f"📍 独立靶区 {locus['locusId']} (定位: {locus['anchorStart']} 附近) - 包含 {len(locus['variants'])} 套微调变体"):
                        for v_idx, cand in enumerate(locus['variants']):
                            role = "🏅 主力优选" if v_idx == 0 else f"🥈 备选 {v_idx}"
                            st.markdown(f"**{role} | 综合得分: {cand['score']:.1f}** (产物长度: {cand['size']} bp | 起始: {cand['start']})")
                            st.code(
                                f"Fwd: {', '.join(cand['fwd'])}  (0错配: {cand['fStats']['p0']}%)\n"
                                f"Prb: {', '.join(cand['probe'])}  (0错配: {cand['pStats']['p0']}%)\n"
                                f"Rev: {', '.join(cand['rev'])}  (0错配: {cand['rStats']['p0']}%)",
                                language="text"
                            )
