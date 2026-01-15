import argparse
import sys
import os

def load_ortholog_map(ortho_file):
    """
    加载直系同源基因映射字典。
    文件格式：GeneA (Species1) \t GeneB (Species2)
    
    :param ortho_file: 直系同源基因文件路径
    :return: 字典 {物种1基因ID: 物种2基因ID}
    """
    ortholog_map = {}
    print(f"[信息] 正在加载直系同源文件: {ortho_file}")
    
    try:
        with open(ortho_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) < 2:
                    # 仅在调试时开启此警告，避免日志过多
                    # print(f"[警告] 第 {line_num} 行格式异常，跳过: {line}")
                    continue
                
                # 建立映射: Lso_Gene -> Lma_Gene
                species1_gene = parts[0]
                species2_gene = parts[1]
                ortholog_map[species1_gene] = species2_gene
                
        print(f"[信息] 直系同源映射加载完成，共 {len(ortholog_map)} 条记录。")
        return ortholog_map

    except FileNotFoundError:
        print(f"[错误] 未找到文件: {ortho_file}")
        sys.exit(1)
    except Exception as e:
        print(f"[错误] 读取文件失败: {e}")
        sys.exit(1)

def process_paralogs(paralog_file, ortholog_map, output_file):
    """
    读取旁系同源文件，结合直系同源映射，生成四元组文件。
    
    输入(Paralog)格式: ID \t GeneA1 \t GeneA2 ...
    输出(Quartet)格式: GeneA1 \t GeneB1 \t GeneA2 \t GeneB2
    
    逻辑:
    1. 读取 Lso vs Lso 的旁系对 (Lso1, Lso2)
    2. 查找 Lso1 对应的 Lma1 (直系)
    3. 查找 Lso2 对应的 Lma2 (直系)
    4. 确保 Lma1 != Lma2 (指向不同的直系基因)
    """
    print(f"[信息] 正在处理旁系同源文件: {paralog_file}")
    
    valid_count = 0
    total_lines = 0
    
    try:
        with open(paralog_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            for line in f_in:
                line = line.strip()
                if not line:
                    continue
                
                total_lines += 1
                parts = line.split()
                
                # 确保至少有3列 (ID, Gene1, Gene2)
                if len(parts) < 3:
                    continue
                
                # 提取旁系基因对 (索引1和2)
                p_gene1 = parts[1]
                p_gene2 = parts[2]
                
                # 核心筛选逻辑：
                # 1. 两个旁系基因必须都在直系同源字典中存在
                if p_gene1 in ortholog_map and p_gene2 in ortholog_map:
                    
                    # 获取对应的直系基因
                    o_gene1 = ortholog_map[p_gene1]
                    o_gene2 = ortholog_map[p_gene2]
                    
                    # 2. 简单的格式校验 (保留原脚本逻辑：需包含 'g' 字符)
                    if 'g' not in o_gene1 or 'g' not in o_gene2:
                        continue
                        
                    # 3. 排除直系基因相同的情况 (防止多对一映射导致的重复)
                    if o_gene1 != o_gene2:
                        # 写入四元组: Lso1, Lma1, Lso2, Lma2
                        f_out.write(f"{p_gene1}\t{o_gene1}\t{p_gene2}\t{o_gene2}\n")
                        valid_count += 1

        print(f"[信息] 处理完成。")
        print(f"[统计] 扫描旁系记录: {total_lines} 条")
        print(f"[统计] 生成有效四元组: {valid_count} 个")
        print(f"[结果] 输出文件已保存至: {output_file}")

    except FileNotFoundError:
        print(f"[错误] 未找到文件: {paralog_file}")
        sys.exit(1)
    except Exception as e:
        print(f"[错误] 处理过程发生异常: {e}")
        sys.exit(1)

def main():
    # 定义命令行参数
    parser = argparse.ArgumentParser(
        description="基于直系同源和旁系同源文件构建基因四元组 (Quartet)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '-o', '--ortho',
        default='Lso_Lma.block.rr.txt.pseu.ortologs',
        help='直系同源基因文件路径 (默认: Lso_Lma.block.rr.txt.pseu.ortologs)'
    )
    parser.add_argument(
        '-p', '--para',
        default='Lso.v.Lso.paralog',
        help='旁系同源基因文件路径 (默认: Lso.v.Lso.paralog)'
    )
    parser.add_argument(
        '-out', '--output',
        default='Lso_Lma.quartet',
        help='输出的四元组文件路径 (默认: Lso_Lma.quartet)'
    )
    
    args = parser.parse_args()
    
    # 步骤 1: 加载映射
    ortho_map = load_ortholog_map(args.ortho)
    
    # 步骤 2: 筛选并输出结果
    process_paralogs(args.para, ortho_map, args.output)

if __name__ == "__main__":
    main()