import sys
import os
import argparse
import subprocess
import random
import math
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

# =========================================================
# 1. 科学计算核心模块 (Ka/Ks 计算)
# =========================================================

# 获取标准遗传密码表 (NCBI Table 1)
STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]

def get_neighbors_score(codon, aa):
    """
    计算一个密码子的同义(Synonymous)和非同义(Nonsynonymous)邻居分数。
    原理：Nei-Gojobori (1986) 算法的简化步骤。
    """
    syn_neighbors = 0
    neighbors = 0
    bases = ['T', 'C', 'A', 'G']
    
    # 对密码子的3个位置分别尝试突变
    for pos in range(3):
        orig_base = codon[pos]
        for b in bases:
            if b == orig_base:
                continue
            
            # 构建突变密码子
            temp_codon_list = list(codon)
            temp_codon_list[pos] = b
            temp_codon = "".join(temp_codon_list)
            
            # 查表获取氨基酸
            # 注意：Biopython的表不包含终止密码子作为键，需捕获异常或预处理
            try:
                temp_aa = STANDARD_TABLE.forward_table.get(temp_codon, '*')
            except KeyError:
                if temp_codon in STANDARD_TABLE.stop_codons:
                    temp_aa = '*'
                else:
                    temp_aa = 'X'

            if temp_aa == '*': continue # 忽略导致终止的突变
            
            neighbors += 1
            if temp_aa == aa:
                syn_neighbors += 1
                
    s_score = syn_neighbors / 3.0
    n_score = (neighbors - syn_neighbors) / 3.0
    return s_score, n_score

def count_sites(seq):
    """
    计算序列中的潜在同义位点(S_sites)和非同义位点(N_sites)。
    """
    S_sites = 0
    N_sites = 0
    length = len(seq)
    
    for i in range(0, length, 3):
        codon = seq[i:i+3]
        if len(codon) < 3 or 'N' in codon or '-' in codon:
            continue
        
        # 获取当前氨基酸
        if codon in STANDARD_TABLE.stop_codons: continue
        aa = STANDARD_TABLE.forward_table.get(codon, 'X')
        if aa == 'X': continue

        s_score, n_score = get_neighbors_score(codon, aa)
        S_sites += s_score
        N_sites += n_score

    return S_sites, N_sites

def calculate_kaks(seq1, seq2):
    """
    计算两序列间的 Ka (非同义替代率) 和 Ks (同义替代率)。
    使用 Jukes-Cantor 校正模型。
    """
    if len(seq1) != len(seq2): return 0, 0, 0, 0
    
    # 1. 计算潜在位点总数
    S_sites1, N_sites1 = count_sites(seq1)
    S_sites2, N_sites2 = count_sites(seq2)
    S_total = (S_sites1 + S_sites2) / 2.0
    N_total = (N_sites1 + N_sites2) / 2.0
    
    if S_total == 0 or N_total == 0: return 0, 0, 0, 0

    sd = 0 # 观察到的同义差异
    nd = 0 # 观察到的非同义差异
    
    # 2. 逐密码子比对
    for i in range(0, len(seq1), 3):
        c1 = seq1[i:i+3]
        c2 = seq2[i:i+3]
        
        if '-' in c1 or '-' in c2 or 'N' in c1 or 'N' in c2: continue
        
        # 跳过终止密码子
        if c1 in STANDARD_TABLE.stop_codons or c2 in STANDARD_TABLE.stop_codons: continue
        
        aa1 = STANDARD_TABLE.forward_table.get(c1, 'X')
        aa2 = STANDARD_TABLE.forward_table.get(c2, 'X')
        
        if aa1 == 'X' or aa2 == 'X': continue
        
        diffs = sum(1 for j in range(3) if c1[j] != c2[j])
        if diffs == 0: continue
        
        # 简化判定：氨基酸改变即为非同义差异
        if aa1 == aa2: sd += 1 
        else: nd += 1 

    # 3. 计算比率 (P-distance)
    Ps = sd / S_total if S_total > 0 else 0
    Pn = nd / N_total if N_total > 0 else 0
    
    # 4. Jukes-Cantor 校正
    def jc_correction(p):
        if p >= 0.75: return 5.0 # 饱和上限
        try:
            return -0.75 * math.log(1 - (4.0/3.0) * p)
        except ValueError:
            return 5.0 # 处理 log 负数情况

    Ks = jc_correction(Ps)
    Ka = jc_correction(Pn)
    
    return Ka, Ks, Pn, Ps

# =========================================================
# 2. 序列处理与比对工具
# =========================================================

def check_dependencies():
    """检查 ClustalW 是否存在"""
    if not shutil.which("clustalw2"):
        print("[错误] 未在系统路径中找到 'clustalw2'。")
        print("       请安装 ClustalW2 (sudo apt install clustalw 或下载二进制包)。")
        sys.exit(1)

def load_fasta(fasta_file):
    """加载FASTA文件到内存"""
    seqs = {}
    if not os.path.exists(fasta_file):
        print(f"[错误] 找不到文件: {fasta_file}")
        sys.exit(1)
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seqs[record.id] = str(record.seq).upper()
    except Exception as e:
        print(f"[错误] 读取 FASTA 文件失败: {e}")
        sys.exit(1)
    return seqs

def run_alignment_workflow(gene_ids, all_seqs, temp_dir):
    """
    执行：提取DNA -> 翻译蛋白 -> ClustalW比对 -> 回译DNA比对
    """
    prot_fasta = os.path.join(temp_dir, "temp_prot.fasta")
    prot_aln = os.path.join(temp_dir, "temp_prot.aln")
    
    prot_records = []
    gene_dna_map = {}
    
    # 1. 翻译
    for uid in gene_ids:
        if uid not in all_seqs:
            return None # 序列缺失
        
        dna = all_seqs[uid]
        try:
            # table=1 (标准), cds=False (允许非完整CDS), to_stop=True (遇终止子停止)
            prot_seq = str(Seq(dna).translate(table=1, to_stop=True))
            if len(prot_seq) < 5: # 忽略极短序列
                return None
            
            prot_records.append(SeqRecord(Seq(prot_seq), id=uid, description=""))
            gene_dna_map[uid] = dna
        except Exception:
            return None

    SeqIO.write(prot_records, prot_fasta, "fasta")
    
    # 2. 调用 ClustalW
    cmd = ["clustalw2", f"-infile={prot_fasta}", f"-outfile={prot_aln}", "-output=FASTA", "-quiet"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return None

    # 3. 回译 (Back-translation)
    if not os.path.exists(prot_aln): return None
    
    dna_aln = {}
    for prot in SeqIO.parse(prot_aln, "fasta"):
        gene_id = prot.id
        aa_seq = str(prot.seq)
        dna_seq = gene_dna_map[gene_id]
        
        codon_aln = []
        dna_idx = 0
        
        for aa in aa_seq:
            if aa == '-':
                codon_aln.append('---')
            else:
                codon_aln.append(dna_seq[dna_idx:dna_idx+3])
                dna_idx += 3
        
        dna_aln[gene_id] = "".join(codon_aln)
        
    return dna_aln

def bootstrap_resample(dna_alignment, length_codons):
    """Bootstrap 重采样"""
    ids = list(dna_alignment.keys())
    new_aln = {i: "" for i in ids}
    cols = list(range(length_codons))
    
    # 有放回随机抽样
    sampled_cols = random.choices(cols, k=length_codons)
    
    for idx in sampled_cols:
        start = idx * 3
        for gene_id in ids:
            new_aln[gene_id] += dna_alignment[gene_id][start:start+3]
            
    return new_aln

# =========================================================
# 3. 主程序
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="基于四元组(Quartet)检测基因转换事件 (Gene Conversion Detection)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-q", "--quartet", required=True, help="基因四元组文件 (由上一步脚本生成)")
    parser.add_argument("-a", "--fasta1", required=True, help="物种1的 CDS 序列文件 (.fasta)")
    parser.add_argument("-b", "--fasta2", required=True, help="物种2的 CDS 序列文件 (.fasta)")
    parser.add_argument("-o", "--output", required=True, help="输出结果文件路径")
    parser.add_argument("--boot", type=int, default=100, help="Bootstrap 重采样次数 (默认: 100)")
    
    args = parser.parse_args()
    
    # 环境检查
    check_dependencies()
    
    # 创建临时目录
    temp_dir = "temp_conversion_work"
    if not os.path.exists(temp_dir): os.makedirs(temp_dir)

    try:
        # 1. 加载数据
        print(f"[信息] 正在加载序列数据...")
        seqs1 = load_fasta(args.fasta1)
        seqs2 = load_fasta(args.fasta2)
        all_seqs = {**seqs1, **seqs2} # 合并字典
        print(f"[信息] 序列加载完成 (物种1序列数: {len(seqs1)}, 物种2序列数: {len(seqs2)})")

        # 2. 准备输出
        with open(args.output, 'w', encoding='utf-8') as out_fh:
            headers = [
                "QuartetID", 
                "Ka_P1", "Ks_P1", "Ka_P2", "Ks_P2", 
                "Ka_O1", "Ks_O1", "Ka_O2", "Ks_O2",
                "Conv_Sp1(Para<Ortho)", "Conv_Sp2(Para<Ortho)", 
                "Boot_Prob_Sp1", "Boot_Prob_Sp2"
            ]
            out_fh.write("\t".join(headers) + "\n")
            
            processed_count = 0
            
            print(f"[信息] 开始分析四元组文件: {args.quartet}")
            
            with open(args.quartet, 'r', encoding='utf-8') as qf:
                for line in qf:
                    line = line.strip()
                    if not line or line.startswith("#"): continue
                    
                    parts = line.split()
                    if len(parts) < 4: continue
                    
                    # 假定格式: Para1 Ortho1 Para2 Ortho2
                    # 对应: Lso1 Lma1 Lso2 Lma2
                    id_p1_a, id_o1, id_p1_b, id_o2 = parts[0], parts[1], parts[2], parts[3]
                    
                    current_ids = [id_p1_a, id_p1_b, id_o1, id_o2]
                    quartet_id = f"{id_p1_a}-{id_p1_b}"
                    
                    # 执行比对流程
                    dna_aln = run_alignment_workflow(current_ids, all_seqs, temp_dir)
                    if not dna_aln:
                        continue
                        
                    # 计算距离
                    # P1: 物种1内的旁系 (Lso1 vs Lso2)
                    # P2: 物种2内的旁系 (Lma1 vs Lma2) - 注意：原脚本逻辑需要两个物种的旁系对比
                    # O1: 直系对1 (Lso1 vs Lma1)
                    # O2: 直系对2 (Lso2 vs Lma2)
                    
                    # 映射关系:
                    # id_p1_a (Lso1)
                    # id_p1_b (Lso2)
                    # id_o1   (Lma1)
                    # id_o2   (Lma2)
                    
                    # 定义对比组
                    pairs = {
                        "P1": (id_p1_a, id_p1_b), # Sp1 Paralog
                        "P2": (id_o1, id_o2),     # Sp2 Paralog (Inferred)
                        "O1": (id_p1_a, id_o1),   # Ortho pair 1
                        "O2": (id_p1_b, id_o2)    # Ortho pair 2
                    }
                    
                    stats = {}
                    valid_calc = True
                    for k, (u1, u2) in pairs.items():
                        res = calculate_kaks(dna_aln[u1], dna_aln[u2])
                        if res[1] == 0 and res[3] == 0: # Ks=0且Ps=0，可能是完全相同或错误
                             # 允许完全相同的情况存在
                             pass
                        stats[k] = res

                    # 提取 Ks 值
                    ks_p1 = stats["P1"][1]
                    ks_p2 = stats["P2"][1]
                    ks_o1 = stats["O1"][1]
                    ks_o2 = stats["O2"][1]
                    
                    # 基因转换判定逻辑 (Ks_Paralog < Ks_Ortholog)
                    # 物种1是否发生转换：Lso1和Lso2比直系更像
                    is_conv_sp1 = (ks_p1 < ks_o1 and ks_p1 < ks_o2)
                    
                    # 物种2是否发生转换：Lma1和Lma2比直系更像
                    is_conv_sp2 = (ks_p2 < ks_o1 and ks_p2 < ks_o2)
                    
                    # Bootstrap 验证
                    boot_sup_sp1 = 0
                    boot_sup_sp2 = 0
                    
                    if is_conv_sp1 or is_conv_sp2:
                        aln_len = len(list(dna_aln.values())[0]) // 3
                        for _ in range(args.boot):
                            res_aln = bootstrap_resample(dna_aln, aln_len)
                            
                            # 快速计算无需全部计算，只需计算必要的
                            b_ks_p1 = calculate_kaks(res_aln[pairs["P1"][0]], res_aln[pairs["P1"][1]])[1]
                            b_ks_p2 = calculate_kaks(res_aln[pairs["P2"][0]], res_aln[pairs["P2"][1]])[1]
                            b_ks_o1 = calculate_kaks(res_aln[pairs["O1"][0]], res_aln[pairs["O1"][1]])[1]
                            b_ks_o2 = calculate_kaks(res_aln[pairs["O2"][0]], res_aln[pairs["O2"][1]])[1]
                            
                            if is_conv_sp1 and (b_ks_p1 < b_ks_o1 and b_ks_p1 < b_ks_o2):
                                boot_sup_sp1 += 1
                            if is_conv_sp2 and (b_ks_p2 < b_ks_o1 and b_ks_p2 < b_ks_o2):
                                boot_sup_sp2 += 1
                                
                    prob_sp1 = boot_sup_sp1 / args.boot if is_conv_sp1 else 0.0
                    prob_sp2 = boot_sup_sp2 / args.boot if is_conv_sp2 else 0.0
                    
                    # 写入行
                    row = [
                        quartet_id,
                        f"{stats['P1'][0]:.4f}", f"{stats['P1'][1]:.4f}",
                        f"{stats['P2'][0]:.4f}", f"{stats['P2'][1]:.4f}",
                        f"{stats['O1'][0]:.4f}", f"{stats['O1'][1]:.4f}",
                        f"{stats['O2'][0]:.4f}", f"{stats['O2'][1]:.4f}",
                        "Y" if is_conv_sp1 else "N",
                        "Y" if is_conv_sp2 else "N",
                        f"{prob_sp1:.2f}",
                        f"{prob_sp2:.2f}"
                    ]
                    out_fh.write("\t".join(row) + "\n")
                    
                    processed_count += 1
                    if processed_count % 10 == 0:
                        print(f"\r[进度] 已处理 {processed_count} 个四元组...", end='', flush=True)

        print(f"\n[完成] 结果已保存至: {args.output}")

    finally:
        # 清理临时目录
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()