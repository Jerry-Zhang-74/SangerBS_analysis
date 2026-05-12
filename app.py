import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import struct
import tempfile
import os

try:
    from snapgene_reader import snapgene_file_to_dict
    HAS_SNAPGENE = True
except ImportError:
    HAS_SNAPGENE = False

st.set_page_config(page_title="Sanger BS测序分析", layout="wide")

PRISM_PASTEL_COLORS = ['#BEBADA', '#FB8072', '#80B1D3', '#8DD3C7', '#FDB462', '#B3DE69', '#FCCDE5', '#FFFFB3', '#D9D9D9', '#BC80BD']

def parse_ab1_custom(file_bytes):
    if file_bytes[:4] != b'ABIF': return None
    root_dir = file_bytes[6:34]
    _, _, _, _, num_elems, data_size, data_offset, _ = struct.unpack('>4sihhiiii', root_dir)
    dir_data = file_bytes[data_offset : data_offset + data_size]
    entries = {}
    for i in range(num_elems):
        entry_bytes = dir_data[i*28 : (i+1)*28]
        if len(entry_bytes) < 28: continue
        e_name, e_tag, e_type, e_size, e_num, e_dsize, e_doffset, _ = struct.unpack('>4sihhiiii', entry_bytes)
        tag_name = f"{e_name.decode('ascii', 'ignore').strip('\x00')}_{e_tag}"
        data_bytes = struct.pack('>i', e_doffset)[:e_dsize] if e_dsize <= 4 else file_bytes[e_doffset : e_doffset + e_dsize]
        entries[tag_name] = {'bytes': data_bytes}
    return entries

def parse_reference(file_obj, filename):
    if filename.lower().endswith(('.fa', '.fasta', '.txt')):
        return "".join([l.strip().upper() for l in file_obj.read().decode('utf-8').splitlines() if l and not l.startswith('>')])
    elif filename.lower().endswith('.dna') and HAS_SNAPGENE:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".dna") as tmp:
            tmp.write(file_obj.read())
            tmp_path = tmp.name
        try:
            return snapgene_file_to_dict(tmp_path)['seq'].upper()
        finally:
            if os.path.exists(tmp_path): os.remove(tmp_path)
    return None

def get_reverse_complement(seq):
    return seq.translate(str.maketrans('ATCGN', 'TAGCN'))[::-1]

def slide_align(ref_3l, ab1_3l):
    best_score = -1
    best_offset = 0
    best_overlap = 0
    len_r = len(ref_3l)
    len_a = len(ab1_3l)
    
    for o in range(-len_r, len_a):
        r_s = max(0, -o)
        r_e = min(len_r, len_a - o)
        overlap = r_e - r_s
        if overlap < 40: continue 
        
        matches = sum(1 for a, b in zip(ref_3l[r_s:r_e], ab1_3l[r_s+o:r_e+o]) if a == b)
        if matches > best_score:
            best_score = matches
            best_offset = o
            best_overlap = overlap
            
    return best_offset, best_score, best_overlap

def analyze_ab1_robust(file_bytes, filename, ref_seq, noise_threshold=5, trim_start=0, trim_end=10, baseline_mode=False):
    try:
        entries = parse_ab1_custom(file_bytes)
        if not entries: return [], None
        
        fwo = entries.get('FWO_1')
        channel_order = "GATC" if not fwo else fwo['bytes'].decode('ascii', 'ignore').replace('\x00', '').strip()
        if len(channel_order) < 4: channel_order = "GATC"
        idx_g, idx_a, idx_t, idx_c = channel_order.find('G'), channel_order.find('A'), channel_order.find('T'), channel_order.find('C')
        
        pbas = entries.get('PBAS_1') or entries.get('PBAS_2')
        ab1_seq = pbas['bytes'].decode('ascii', 'ignore').replace('\x00', '')
        
        ploc_entry = entries.get('PLOC_2') or entries.get('PLOC_1')
        ploc = struct.unpack(f'>{len(ploc_entry["bytes"])//2}h', ploc_entry["bytes"])
        
        def get_trace(idx):
            for i in [9, 1]:
                tag = f"DATA_{i + idx}"
                if tag in entries: return struct.unpack(f'>{len(entries[tag]["bytes"])//2}h', entries[tag]['bytes'])
            return None

        trace_g, trace_a, trace_t, trace_c = get_trace(idx_g), get_trace(idx_a), get_trace(idx_t), get_trace(idx_c)
        if not all([trace_g, trace_a, trace_t, trace_c]): return [], None

        fwd_ref_3l = ref_seq.replace('C', 'T')
        fwd_ab1_3l = ab1_seq.replace('C', 'T')
        fwd_offset, fwd_score, fwd_overlap = slide_align(fwd_ref_3l, fwd_ab1_3l)
        
        ref_rc = get_reverse_complement(ref_seq)
        rev_ref_3l = ref_rc.replace('G', 'A')
        rev_ab1_3l = ab1_seq.replace('G', 'A')
        rev_offset, rev_score, rev_overlap = slide_align(rev_ref_3l, rev_ab1_3l)
        
        is_reverse = rev_score > fwd_score
        best_offset = rev_offset if is_reverse else fwd_offset
        best_score = rev_score if is_reverse else fwd_score
        best_overlap = rev_overlap if is_reverse else fwd_overlap
        
        active_ref = ref_rc if is_reverse else ref_seq
        active_ref_3l = rev_ref_3l if is_reverse else fwd_ref_3l
        active_ab1_3l = rev_ab1_3l if is_reverse else fwd_ab1_3l
        
        identity = (best_score / best_overlap * 100) if best_overlap > 0 else 0
        if identity < 40:
            return [], {"Sample": filename, "QC Status": f"比对失败 (同源率 {identity:.1f}%)"}

        total_non_cpg, unconverted = 0, 0
        unconverted_absolute = 0
        for ref_pos in range(len(active_ref)):
            ab1_pos = ref_pos + best_offset
            if 0 <= ab1_pos < len(ab1_seq) and trim_start <= ab1_pos <= len(ab1_seq) - trim_end:
                if is_reverse:
                    if active_ref[ref_pos] == 'G' and (ref_pos == 0 or active_ref[ref_pos-1] != 'C'):
                        total_non_cpg += 1
                        if ab1_seq[ab1_pos] == 'G': unconverted += 1
                    if ab1_seq[ab1_pos] == 'G' and (ab1_pos == 0 or ab1_seq[ab1_pos-1] != 'C'):
                        unconverted_absolute += 1
                else:
                    if active_ref[ref_pos] == 'C' and (ref_pos == len(active_ref)-1 or active_ref[ref_pos+1] != 'G'):
                        total_non_cpg += 1
                        if ab1_seq[ab1_pos] == 'C': unconverted += 1
                    if ab1_seq[ab1_pos] == 'C' and (ab1_pos == len(ab1_seq)-1 or ab1_seq[ab1_pos+1] != 'G'):
                        unconverted_absolute += 1

        dir_str = "反向" if is_reverse else "正向"
        if total_non_cpg > 0:
            qc_status = f"{dir_str} (同源:{identity:.1f}%) | 转化率 {100.0 - (unconverted / total_non_cpg * 100.0):.1f}%"
        else:
            qc_status = f"{dir_str} (同源:{identity:.1f}%) | 参考已转换(残留未转化C: {unconverted_absolute}个)"

        original_ref_cpgs = [i for i in range(len(ref_seq)-1) if ref_seq[i:i+2] == 'CG']
        results = []

        for i, orig_ref_pos in enumerate(original_ref_cpgs):
            meth = np.nan
            active_ref_pos = len(ref_seq) - 1 - orig_ref_pos if is_reverse else orig_ref_pos
            expected_ab1_pos = active_ref_pos + best_offset 
            
            c_start = max(0, active_ref_pos - 15)
            c_end = min(len(active_ref_3l), active_ref_pos + 17)
            ctx = active_ref_3l[c_start:c_end]
            target_idx_in_ctx = active_ref_pos - c_start
            
            s_start = max(0, expected_ab1_pos - 30)
            s_end = min(len(ab1_seq), expected_ab1_pos + 30)
            
            best_ctx_score = -1
            best_ab1_pos = expected_ab1_pos
            
            if len(ctx) >= 10 and s_start < s_end - len(ctx):
                for j in range(s_start, s_end - len(ctx) + 1):
                    score = sum(1 for a, b in zip(ctx, active_ab1_3l[j:j+len(ctx)]) if a == b)
                    if score > best_ctx_score:
                        best_ctx_score = score
                        best_ab1_pos = j + target_idx_in_ctx
                        
            final_ab1_pos = best_ab1_pos if best_ctx_score >= len(ctx) * 0.65 else expected_ab1_pos

            if 0 <= final_ab1_pos < len(ab1_seq) and trim_start <= final_ab1_pos <= len(ab1_seq) - trim_end:
                if final_ab1_pos < len(ploc):
                    loc = ploc[final_ab1_pos]
                    
                    t_meth = trace_g if is_reverse else trace_c
                    t_unmeth = trace_a if is_reverse else trace_t

                    win_s, win_e = max(0, loc-8), min(len(t_meth), loc+9)
                    if win_e > win_s:
                        sum_trace = [t_meth[k] + t_unmeth[k] for k in range(win_s, win_e)]
                        center_idx = win_s + np.argmax(sum_trace)

                        c_s, c_e = max(0, center_idx-2), min(len(t_meth), center_idx+3)
                        t_s, t_e = max(0, center_idx-2), min(len(t_unmeth), center_idx+3)
                        
                        h_m_raw = max(t_meth[c_s:c_e]) if c_e > c_s else 0
                        h_u_raw = max(t_unmeth[t_s:t_e]) if t_e > t_s else 0

                        if baseline_mode:
                            bg_m = min(t_meth[max(0, center_idx-10):min(len(t_meth), center_idx+11)])
                            bg_u = min(t_unmeth[max(0, center_idx-10):min(len(t_unmeth), center_idx+11)])
                            h_m = max(0, h_m_raw - bg_m)
                            h_u = max(0, h_u_raw - bg_u)
                        else:
                            h_m, h_u = h_m_raw, h_u_raw

                        if (h_m + h_u) >= noise_threshold:
                            meth = (h_m / (h_m + h_u)) * 100.0

            results.append({
                "Sample": filename,
                "CpG_Site": f"Site_{i+1}",
                "Ref_Pos": orig_ref_pos + 1,
                "Methylation (%)": float(round(meth, 2)) if not np.isnan(meth) else np.nan
            })

        return results, {"Sample": filename, "QC Status": qc_status}
    except Exception as e:
        st.error(f"解析出错: {e}")
        return [], None

st.title("Sanger BS测序分析工具")

with st.sidebar:
    st.header("1. 参考序列输入")
    ref_file = st.file_uploader("上传参考序列 (.dna/.fasta)", type=['dna', 'fasta', 'fa', 'txt'])
    
    st.header("2. 分析参数设置")
    trim_start = st.number_input("忽略前端碱基数 (bp)", value=0, step=5)
    trim_end = st.number_input("忽略尾端碱基数 (bp)", value=10, step=5)
    noise_threshold = st.number_input("信号总和最低阈值", value=5, step=5)
    baseline_mode = st.checkbox("扣除局部基线", value=False)
    
    st.header("3. 测序文件导入")
    uploaded_files = st.file_uploader("上传 .ab1 测序文件集", type=['ab1'], accept_multiple_files=True)

if uploaded_files:
    if not ref_file:
        st.warning("请先上传参考序列以进行比对和位点识别。")
    else:
        ref_seq = parse_reference(ref_file, ref_file.name)
        st.info(f"参考序列已加载。共识别出 {ref_seq.count('CG')} 个 CpG 位点。")
        
        all_results, all_qc, seen_names = [], [], {}
        with st.spinner("正在执行序列比对与信号提取，请稍候..."):
            for file in uploaded_files:
                base = file.name
                u_name = f"{base} ({seen_names[base]})" if base in seen_names else base
                seen_names[base] = seen_names.get(base, 0) + 1
                    
                res, qc = analyze_ab1_robust(file.read(), u_name, ref_seq, noise_threshold, trim_start, trim_end, baseline_mode)
                if res: all_results.extend(res)
                if qc: all_qc.append(qc)
                
        if all_qc:
            st.subheader("实验质控 (QC)")
            st.dataframe(pd.DataFrame(all_qc), width='stretch')
            st.divider()

        if all_results:
            df = pd.DataFrame(all_results)
            
            with st.expander("点击展开：查看原始明细数据"):
                st.dataframe(df, width='stretch')
            
            try:
                df['Methylation (%)'] = pd.to_numeric(df['Methylation (%)'], errors='coerce')
                df['Site_Num'] = df['CpG_Site'].str.extract(r'(\d+)').astype(int)
                
                df = df.drop_duplicates(subset=['Sample', 'CpG_Site', 'Ref_Pos'])
                pivot = df.set_index(['CpG_Site', 'Ref_Pos', 'Sample'])['Methylation (%)'].unstack('Sample')
                
                site_map = df.drop_duplicates('CpG_Site').set_index('CpG_Site')['Site_Num']
                pivot['Site_Num'] = pivot.index.get_level_values('CpG_Site').map(site_map)
                pivot = pivot.sort_values('Site_Num').drop(columns=['Site_Num'])
                
                st.subheader("甲基化分析结果矩阵")
                st.dataframe(pivot, width='stretch')
                
                col1, col2 = st.columns(2)
                col1.download_button("下载分析矩阵 (CSV)", data=pivot.to_csv().encode('utf-8-sig'), file_name="Methylation_Matrix.csv", mime="text/csv", width='stretch')
                col2.download_button("下载明细数据 (CSV)", data=df.to_csv(index=False).encode('utf-8-sig'), file_name="Methylation_Details.csv", mime="text/csv", width='stretch')
            except Exception as pivot_e:
                st.error(f"矩阵生成发生错误: {pivot_e}")
            
            st.subheader("甲基化水平折线图")
            plot_df = df.dropna(subset=['Methylation (%)'])
            if not plot_df.empty:
                fig = go.Figure()
                for i, sample in enumerate(plot_df['Sample'].unique()):
                    s_data = plot_df[plot_df['Sample'] == sample].sort_values('Ref_Pos')
                    fig.add_trace(go.Scatter(
                        x=s_data['Ref_Pos'], y=s_data['Methylation (%)'],
                        mode='lines+markers', name=sample,
                        line=dict(width=2.5, color=PRISM_PASTEL_COLORS[i % 10]),
                        marker=dict(size=9, line=dict(width=0)),
                        customdata=s_data['CpG_Site'],
                        hovertemplate='<b>%{customdata}</b><br>Pos: %{x} bp<br>Methylation: %{y:.1f}%<extra></extra>'
                    ))
                fig.update_layout(height=500, plot_bgcolor='rgba(0,0,0,0)', hovermode='x unified',
                                  yaxis=dict(title="Methylation (%)", range=[-5, 105], showline=True, linewidth=2, linecolor='black'),
                                  xaxis=dict(title="Reference Sequence Position (bp)", showline=True, linewidth=2, linecolor='black', rangeslider=dict(visible=True, thickness=0.15)),
                                  legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
                st.plotly_chart(fig, width='stretch')
