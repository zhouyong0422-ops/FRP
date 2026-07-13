import streamlit as st
import csv
import re
import io
import math
from collections import Counter
from datetime import datetime

# ==============================================================================
# 页面基础与 UI 配置
# ==============================================================================
st.set_page_config(page_title="V02 原料齐套 - 工业级全景长度竞争寻优系统 V16", layout="wide")

# ==============================================================================
# 底层核心计算函数 (无损复刻 V16 HTML 版逻辑)
# ==============================================================================
def read_fasta_from_string(fasta_string):
    """读取并解析 FASTA 序列流"""
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
    """计算经验 Tm 值"""
    g = seq.count('G')
    c = seq.count('C')
    return round(64.9 + 41 * (g + c - 16.4) / len(seq), 1)

def calc_gc(seq):
    """计算 GC 含量 (%)"""
    g = seq.count('G')
    c = seq.count('C')
    return round(((g + c) / len(seq)) * 100, 1)

def reverse_complement(seq):
    """反向互补"""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(mapping.get(b, b) for b in reversed(seq))

def has_secondary_structure_risk(seq):
    """排查二级结构与发夹自互补风险"""
    length = len(seq)
    if length < 12:
        return False
    for i in range(4, 6):
        head = seq[:i]
        tail_comp = reverse_complement(seq[-i:])
        if head == tail_comp:
            return True
    return False

def has_3prime_dimer_risk(seq1, seq2):
    """3' 端交叉二聚体排查：末端 4bp 绝对不能互补"""
    end1 = seq1[-4:]
    end2_comp = reverse_complement(seq2[-4:])
    return end1 == end2_comp

# 1. 绝对硬底线质控（触犯一票否决）
def is_hard_valid_oligo(seq, is_probe=False):
    if 'N' in seq or re.search(r'[^ATGC]', seq):
        return False
    gc = float(calc_gc(seq))
    if gc < 20 or gc > 80:
        return False
    if is_probe and seq[0] == 'G':
        return False
    return True

# 2. 柔性优先准则（触犯则执行降维惩罚打分）
def calc_soft_penalties(seq, is_probe=False):
    penalty = 0.0
    gc = float(calc_gc(seq))
    
    if re.search(r'([ATGC])\1{3,}', seq):
        penalty += 10
    if has_secondary_structure_risk(seq):
        penalty += 15

    if is_probe:
        if gc < 30: penalty += (30 - gc) * 2
        if gc > 65: penalty += (gc - 65) * 2
        g_count = seq.count('G')
        c_count = seq.count('C')
        if g_count >= c_count:
            penalty += 8
    else:
        if gc < 40: penalty += (40 - gc) * 2
        if gc > 60: penalty += (gc - 60) * 2
        
        end5 = seq[-5:]
        end_gc = end5.count('G') + end5.count('C')
        if end_gc < 1 or end_gc > 2:
            penalty += 6
            
        if seq[-1] == 'T':
            penalty += 12
        if re.search(r'GG$|CC$', seq):
            penalty += 8
            
    return penalty

def get_top_variants(start_index, length, sequences_array, max_variants=2):
    """90% 覆盖率硬熔断及变体提取"""
    counts = Counter()
    total_valid = 0
    total_seq = len(sequences_array)

    for seq in sequences_array:
        sub = seq[start_index : start_index + length]
        if '-' in sub or 'N' in sub:
            continue
        counts[sub] += 1
        total_valid += 1

    if total_seq == 0 or total_valid / total_seq < 0.90:
        return []
    if not counts:
        return []

    sorted_counts = counts.most_common()
    variants = [sorted_counts[0][0]]
    coverage = sorted_counts[0][1] / total_valid

    # 仅当第一大序列覆盖率偏低，且第二大序列占比超 4% 时，才允许提取变体
    if coverage < 0.97 and len(sorted_counts) > 1 and max_variants > 1:
        second_cov = sorted_counts[1][1] / total_valid
        if second_cov > 0.04:
            variants.append(sorted_counts[1][0])
            
    return variants

# 【核心重大重构】：废除 break 锁！对所有合法长度进行竞争打分
def get_all_valid_variants(start_idx, sequences, is_probe=False):
    valid_list = []
    target_lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30] if is_probe else [18, 19, 20, 21, 22, 23, 24, 25]
    
    for length in target_lengths:
        if start_idx + length > len(sequences[0]):
            continue
        raw_variants = get_top_variants(start_idx, length, sequences, 1 if is_probe else 2)
        if not raw_variants:
            continue

        final_variants = []
        all_passed = True
        
        for v in raw_variants:
            if not is_hard_valid_oligo(v, is_probe):
                all_passed = false if 'false' in locals() else False # Python boolean
                break
            final_variants.append(v)
            
        if all_passed:
            soft_pen = sum(calc_soft_penalties(v, is_probe) for v in final_variants)
            avg_tm = sum(float(calc_tm(v)) for v in final_variants) / len(final_variants)
            tm_bonus = 0
            if is_probe:
                if 62 <= avg_tm <= 68: tm_bonus = 15
            else:
                if 55 <= avg_tm <= 60: tm_bonus = 15

            valid_list.append({
                'length': length,
                'variants': final_variants,
                'score': tm_bonus - soft_pen
            })
            
    valid_list.sort(key=lambda x: x['score'], reverse=True)
    return valid_list

def get_reverse_valid_variants(start_idx, sequences):
    valid_list = []
    target_lengths = [18, 19, 20, 21, 22, 23, 24, 25]
    for length in target_lengths:
        if start_idx + length > len(sequences[0]):
            continue
        raw_variants = get_top_variants(start_idx, length, sequences, 2)
        if not raw_variants:
            continue

        final_variants = []
        all_passed = True
        for rv in raw_variants:
            comp = reverse_complement(rv)
            if not is_hard_valid_oligo(comp, False):
                all_passed = False
                break
            final_variants.append(comp)
            
        if all_passed:
            soft_pen = sum(calc_soft_penalties(v, False) for v in final_variants)
            avg_tm = sum(float(calc_tm(v)) for v in final_variants) / len(final_variants)
            tm_bonus = 15 if (55 <= avg_tm <= 60) else 0
            valid_list.append({
                'length': length,
                'variants': final_variants,
                'rawVariants': raw_variants,
                'score': tm_bonus - soft_pen
            })
            
    valid_list.sort(key=lambda x: x['score'], reverse=True)
    return valid_list

def calculate_mix_mismatch(variants_array, start_index, sequences_array):
    total_seq = len(sequences_array)
    stats = {'m0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': total_seq}
    seq_len = len(variants_array[0])

    for seq in sequences_array:
        lib_seq_snippet = seq[start_index : start_index + seq_len]
        if '-' in lib_seq_snippet or 'N' in lib_seq_snippet:
            stats['m3p'] += 1
            continue
        
        best_mismatches = seq_len
        for target_seq in variants_array:
            mismatches = sum(1 for a, b in zip(target_seq, lib_seq_snippet) if a != b)
            if mismatches < best_mismatches:
                best_mismatches = mismatches
                
        if best_mismatches == 0: stats['m0'] += 1
        elif best_mismatches == 1: stats['m1'] += 1
        elif best_mismatches == 2: stats['m2'] += 1
        else: stats['m3p'] += 1

    if stats['total'] == 0:
        return {'p0': '0.0', 'p1': '0.0', 'p2': '0.0', 'p3': '0.0', 'm0': 0, 'm1': 0, 'm2': 0, 'm3p': 0, 'total': 0}
        
    return {
        'p0': round((stats['m0'] / stats['total']) * 100, 1),
        'p1': round((stats['m1'] / stats['total']) * 100, 1),
        'p2': round((stats['m2'] / stats['total']) * 100, 1),
        'p3': round((stats['m3p'] / stats['total']) * 100, 1),
        'm0': stats['m0'], 'm1': stats['m1'], 'm2': stats['m2'], 'm3p': stats['m3p'], 'total': stats['total']
    }

def build_csv_string(global_loci_groups):
    """构建内存 CSV 下载流"""
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["靶区归属", "变体角色", "综合得分", "寡核苷酸类型", "序列 (5'->3')", "长度 (bp)", "Tm (°C)", "GC (%)", "完美匹配(0)", "错配1碱基(1)", "错配2碱基(2)", "错配3碱基(≥3)", "预期产物长度 (bp)", "精确起始坐标"])
    
    for locus in global_loci_groups:
        for v_idx, cand in enumerate(locus['variants']):
            locus_name = f"靶区_{locus['locusId']}"
            role = "主力优选" if v_idx == 0 else f"备选_{v_idx}"
            score = round(cand['score'], 1)
            size = cand['size']
            start = cand['start']
            
            for i, seq in enumerate(cand['fwd']):
                type_str = f"Forward_{i+1}" if len(cand['fwd']) > 1 else "Forward"
                st = cand['fStats']
                writer.writerow([locus_name, role, score, type_str, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])
            
            for i, seq in enumerate(cand['probe']):
                type_str = f"Probe_{i+1}" if len(cand['probe']) > 1 else "Probe"
                st = cand['pStats']
                writer.writerow([locus_name, role, score, type_str, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])

            for i, seq in enumerate(cand['rev']):
                type_str = f"Reverse_{i+1}" if len(cand['rev']) > 1 else "Reverse"
                st = cand['rStats']
                writer.writerow([locus_name, role, score, type_str, seq, len(seq), calc_tm(seq), calc_gc(seq), f"{st['p0']}% ({st['m0']}/{st['total']})", f"{st['p1']}% ({st['m1']}/{st['total']})", f"{st['p2']}% ({st['m2']}/{st['total']})", f"{st['p3']}% ({st['m3p']}/{st['total']})", size, start])
    return output.getvalue()

# ==============================================================================
# Streamlit Web 界面与渲染逻辑
# ==============================================================================
st.title("🧬 自动化引物探针柔性优先寻优系统")
st.markdown("### <span style='background-color:#8e44ad; color:white; padding:4px 12px; border-radius:12px; font-size:14px;'>全量长度竞争 V16</span>", unsafe_allow_html=True)
st.caption("已废除贪心锁定 Bug！引物 18-25nt、探针 18-30nt 全量长度并发竞争，择优录取！死守4大绝对底线；取消短产物偏好；重罚混合套数（极力保障单套纯净输出）。")

uploaded_file = st.file_uploader("📂 导入对齐后的 FASTA 序列库 (支持 fasta, fas, txt, aln)", type=['fasta', 'fas', 'txt', 'aln'])

if uploaded_file is not None:
    fasta_string = uploaded_file.getvalue().decode("utf-8")
    sequences = read_fasta_from_string(fasta_string)
    
    if len(sequences) < 2:
        st.error("❌ 序列读取失败，或文件包含的序列少于 2 条，请检查文件格式。")
    else:
        seq_len = len(sequences[0])
        st.success(f"✅ 成功读取 {len(sequences)} 条对齐序列，对齐总长度：{seq_len} bp。")
        
        # 渲染香农熵变异强度图
        with st.expander("📊 查看靶标序列群变异强度扫描 (香农熵)", expanded=False):
            entropies = []
            for i in range(seq_len):
                column = {}
                total = 0
                for j in range(len(sequences)):
                    base = sequences[j][i]
                    if base and base != '-':
                        column[base] = column.get(base, 0) + 1
                        total += 1
                entropy = 0.0
                for base in column:
                    p = column[base] / total
                    entropy -= p * math.log2(p)
                entropies.append(entropy)
                
            st.bar_chart(entropies)

        if st.button("⚙️ 启动寻优：基于全区长度竞争出具靶区 DOE", type="primary", use_container_width=True):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            all_candidates = []
            min_gap = 1
            max_gap = 50
            total_steps = seq_len - 160
            
            status_text.info("⏳ 正在进行全量长度矩阵竞争与热力学加权计算，请耐心等待...")
            
            for i in range(total_steps):
                if i % 20 == 0 or i == total_steps - 1:
                    progress_bar.progress(int((i / total_steps) * 100))
                    status_text.text(f"🔍 正在全量并发扫描位点 {i} / {total_steps} ...")
                    
                f_obj_list = get_all_valid_variants(i, sequences, False)
                if not f_obj_list: continue

                # 【核心并发点】：取出当前位置得分最高的 Top 2 最优前向引物长度
                for f_obj in f_obj_list[:2]:
                    f_variants = f_obj['variants']
                    f_len = f_obj['length']

                    for gap1 in range(min_gap, max_gap + 1):
                        p_start = i + f_len + gap1
                        p_obj_list = get_all_valid_variants(p_start, sequences, True)
                        if not p_obj_list: continue

                        # 取出当前位置得分最高的 Top 2 最优探针长度
                        for p_obj in p_obj_list[:2]:
                            p_variants = p_obj['variants']
                            p_len = p_obj['length']

                            for gap2 in range(min_gap, max_gap + 1):
                                r_start = p_start + p_len + gap2
                                r_obj_list = get_reverse_valid_variants(r_start, sequences)
                                if not r_obj_list: continue

                                # 取出当前位置得分最高的 Top 2 最优反向引物长度
                                for r_obj in r_obj_list[:2]:
                                    r_variants_raw = r_obj['rawVariants']
                                    r_variants = r_obj['variants']
                                    r_len = r_obj['length']

                                    amplicon_size = r_start + r_len - i
                                    if amplicon_size < 70 or amplicon_size > 160:
                                        continue

                                    min_f_tm = min(float(calc_tm(v)) for v in f_variants)
                                    max_f_tm = max(float(calc_tm(v)) for v in f_variants)
                                    min_r_tm = min(float(calc_tm(v)) for v in r_variants)
                                    max_r_tm = max(float(calc_tm(v)) for v in r_variants)
                                    min_p_tm = min(float(calc_tm(v)) for v in p_variants)
                                    
                                    primer_max_tm = max(max_f_tm, max_r_tm)
                                    primer_min_tm = min(min_f_tm, min_r_tm)
                                    
                                    soft_penalty = 0.0
                                    for f in f_variants: soft_penalty += calc_soft_penalties(f, False)
                                    for p in p_variants: soft_penalty += calc_soft_penalties(p, True)
                                    for r in r_variants: soft_penalty += calc_soft_penalties(r, False)

                                    dimer_risk = any(has_3prime_dimer_risk(f, r) for f in f_variants for r in r_variants)
                                    if dimer_risk: soft_penalty += 20

                                    primer_tm_diff = abs(primer_max_tm - primer_min_tm)
                                    if primer_tm_diff > 2.0: soft_penalty += (primer_tm_diff - 2.0) * 5

                                    if min_p_tm < primer_max_tm + 5.0:
                                        soft_penalty += (primer_max_tm + 5.0 - min_p_tm) * 6

                                    if gap1 > 5: soft_penalty += (gap1 - 5) * 2
                                    if gap2 > 5: soft_penalty += (gap2 - 5) * 2

                                    f_stats = calculate_mix_mismatch(f_variants, i, sequences)
                                    p_stats = calculate_mix_mismatch(p_variants, p_start, sequences)
                                    r_stats = calculate_mix_mismatch(r_variants_raw, r_start, sequences)

                                    f_p0 = float(f_stats['p0'])
                                    p_p0 = float(p_stats['p0'])
                                    r_p0 = float(r_stats['p0'])

                                    base_score = f_p0 + (p_p0 * 3) + r_p0
                                    probe_bonus = 50 if (p_p0 >= 99.0) else 0
                                    probe_penalty = (98.0 - p_p0) * 10 if (p_p0 < 98.0) else 0

                                    mix_f = -35 if (len(f_variants) > 1) else 0
                                    mix_r = -35 if (len(r_variants) > 1) else 0
                                    
                                    total_score = base_score + probe_bonus - probe_penalty + mix_f + mix_r - soft_penalty

                                    all_candidates.append({ 
                                        'fwd': f_variants, 'rev': r_variants, 'probe': p_variants, 
                                        'fStats': f_stats, 'pStats': p_stats, 'rStats': r_stats,
                                        'size': amplicon_size, 'start': i, 'score': total_score,
                                        'details': {'base': base_score, 'pBonus': probe_bonus, 'pPenalty': -probe_penalty, 'mixF': mix_f, 'mixR': mix_r, 'softPen': -soft_penalty}
                                    })

            progress_bar.empty()
            status_text.empty()
            
            # 空间非极大值抑制 (NMS) 聚类
            all_candidates.sort(key=lambda x: x['score'], reverse=True)
            global_loci_groups = []
            locus_window = 50 

            for cand in all_candidates:
                found_locus = False
                for locus in global_loci_groups:
                    if abs(cand['start'] - locus['anchorStart']) <= locus_window:
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

            # 结果生成与展示
            if not global_loci_groups:
                st.error("⚠️ **体系设计失败**：在该序列库中未能找到产物 70-160bp、探针不以 G 开头且覆盖度达标的黄金区。")
            else:
                st.success(f"🎉 **寻优完成！** 共提炼出 {len(global_loci_groups)} 个独立黄金靶区。")
                
                csv_str = build_csv_string(global_loci_groups)
                st.download_button(
                    label="📥 一键导出完整 DOE 清单 (Excel CSV)",
                    data=csv_str.encode('utf-8-sig'),
                    file_name=f"V02_全量长度竞争优选_DOE清单_{datetime.now().strftime('%Y-%m-%d')}.csv",
                    mime="text/csv",
                    type="primary"
                )
                
                st.markdown("---")
                st.subheader("🎯 独立黄金靶区报告 (Top 靶区展示)")

                for locus in global_loci_groups:
                    with st.container():
                        st.markdown(f"#### 📍 独立黄金靶区 {locus['locusId']} <span style='font-size:14px; color:#7f8c8d; font-weight:normal;'> (参考起始坐标: {locus['anchorStart']})</span>", unsafe_allow_html=True)
                        
                        for v_idx, cand in enumerate(locus['variants']):
                            is_primary = (v_idx == 0)
                            role_text = "主力优选" if is_primary else f"微调备选 {v_idx}"
                            role_color = "#27ae60" if is_primary else "#7f8c8d"
                            card_border = "#27ae60" if is_primary else "#bdc3c7"
                            card_bg = "#f9fbfd" if is_primary else "#ffffff"
                            
                            with st.expander(f"🏅 [{role_text}] 综合得分: {cand['score']:.1f} 分 | 精确定位: {cand['start']} | 产物长度: {cand['size']} bp", expanded=is_primary):
                                
                                # 序列基础与统计函数渲染器
                                def render_oligo_block(title, variants_list, stats, is_probe=False):
                                    border_style = "border-left: 4px solid #8e44ad; background: #fdfafb;" if is_probe else "border: 1px solid #ecf0f1; background: #fff;"
                                    html_content = f"<div style='margin-bottom: 10px; padding: 10px; border-radius: 6px; {border_style}'>"
                                    for idx, v_seq in enumerate(variants_list):
                                        label = f"{title} {idx+1}" if len(variants_list) > 1 else title
                                        html_content += f"<div style='display:flex; justify-content:space-between; font-family:monospace; margin-bottom:4px;'><span style='font-weight:bold; color:#34495e; width:90px;'>{label}:</span><span style='color:#d35400; font-weight:bold; flex-grow:1;'>5'- {v_seq} -3'</span><span style='color:#7f8c8d; font-size:12px;'>Len: {len(v_seq)}bp | Tm: {calc_tm(v_seq)}°C | GC: {calc_gc(v_seq)}%</span></div>"
                                    
                                    is_mix_tag = "<span style='background:#9b59b6; color:white; padding:2px 8px; border-radius:8px; font-size:11px; margin-left:8px;'>混合套数扣分(-35)</span>" if len(variants_list) > 1 else ""
                                    html_content += f"<div style='display:flex; flex-wrap:wrap; gap:8px; font-size:12px; margin-top:8px; padding-top:6px; border-top:1px dashed #ecf0f1;'><span style='background:#d4edda; color:#155724; padding:2px 8px; border-radius:10px; border:1px solid #c3e6cb;'>完全匹配(0): {stats['p0']}%, {stats['m0']}/{stats['total']}</span><span style='background:#fff3cd; color:#856404; padding:2px 8px; border-radius:10px; border:1px solid #ffeeba;'>错配1碱基(1): {stats['p1']}%, {stats['m1']}/{stats['total']}</span><span style='background:#ffeeba; color:#856404; padding:2px 8px; border-radius:10px; border:1px solid #ffdf7e;'>错配2碱基(2): {stats['p2']}%, {stats['m2']}/{stats['total']}</span><span style='background:#f8d7da; color:#721c24; padding:2px 8px; border-radius:10px; border:1px solid #f5c6cb;'>错配3碱基(≥3): {stats['p3']}%, {stats['m3p']}/{stats['total']}</span>{is_mix_tag}</div></div>"
                                    st.markdown(html_content, unsafe_allow_html=True)

                                render_oligo_block("Forward", cand['fwd'], cand['fStats'])
                                render_oligo_block("Probe", cand['probe'], cand['pStats'], is_probe=True)
                                render_oligo_block("Reverse", cand['rev'], cand['rStats'])
                                
                                # 明细展示
                                st.markdown("---")
                                d = cand['details']
                                st.markdown(f"""
                                **🔍 柔性评分核算明细：**
                                * **基础匹配分** *(探针3倍权重)*: `+{d['base']:.1f}`
                                * **探针卓越奖励** *(完美率 ≥99%)*: `+{d['pBonus']:.1f}`
                                * **探针错配惩罚** *(低于98%十倍扣除)*: `{d['pPenalty']:.1f}`
                                * **F/R 混合套数重罚** *(极力优选单套)*: `{d['mixF'] + d['mixR']}`
                                * **柔性偏离总扣分** *(GC/温差/3'末位等)*: `{d['softPen']:.1f}`
                                """)
                        st.write("")
