import glob
import os
import argparse
import re
import sys

# -----------------------------------------------------------------------------
# 全局常量与预编译正则
# -----------------------------------------------------------------------------

# 预编译正则表达式以提高解析性能
# 逻辑：匹配开头的一个或多个字母，捕获紧接着的数字，直到遇到字符 'g'
# 例如：对于 "Lso01g12345"，匹配 "Lso01g"，捕获组1为 "01"
GENE_ID_PATTERN = re.compile(r'^[a-zA-Z]+(\d+)g')

# -----------------------------------------------------------------------------
# 辅助函数
# -----------------------------------------------------------------------------

def parse_chromosome(gene_id):
    """
    使用正则解析基因ID中的染色体编号。
    :param gene_id: 基因ID字符串 (如 Lso01g0001)
    :return: 整数类型的染色体号 (如 1)。若解析失败返回 0。
    """
    match = GENE_ID_PATTERN.match(gene_id)
    if match:
        return int(match.group(1)) # 自动处理前导0，例如 int('01') -> 1
    return 0

def parse_target_chroms(chrom_arg):
    """
    解析命令行传入的染色体对参数。
    :param chrom_arg: 字符串 (如 "1,1;11,11")
    :return: 包含元组的集合 (如 {(1,1), (11,11)})
    """
    target_pairs = set()
    try:
        # 按分号分割多组
        pairs = chrom_arg.split(';')
        for pair in pairs:
            if not pair.strip():
                continue
            # 按逗号分割每组的两个染色体
            c1, c2 = pair.split(',')
            target_pairs.add((int(c1.strip()), int(c2.strip())))
    except ValueError:
        print(f"[错误] 参数格式错误: '{chrom_arg}'。正确格式示例: '1,1;11,11'")
        sys.exit(1)
    return target_pairs

# -----------------------------------------------------------------------------
# 核心处理逻辑
# -----------------------------------------------------------------------------

def process_block_files(input_pattern, output_suffix, target_pairs):
    """
    遍历文件并筛选符合条件的基因对。
    """
    # 获取所有匹配的文件
    files = glob.glob(input_pattern)
    if not files:
        print(f"[警告] 未找到匹配模式 '{input_pattern}' 的文件，请检查路径。")
        return

    print(f"目标染色体对: {target_pairs}")
    print(f"待处理文件数: {len(files)}")

    for in_file in files:
        # 构造输出文件名：原文件名 + 后缀
        out_file = f"{in_file}{output_suffix}"
        
        filtered_count = 0
        line_processed = 0
        
        print(f"--------------------------------------------------")
        print(f"正在处理: {in_file}")

        try:
            with open(in_file, 'r', encoding='utf-8') as f_in, \
                 open(out_file, 'w', encoding='utf-8') as f_out:
                
                for line in f_in:
                    line = line.strip()
                    # 快速跳过空行或非基因数据行（假设基因行以'L'开头）
                    # 根据你的文件样本，block header 是 "the 1th path..." 或 "++++..."
                    if not line or not line.startswith('L'):
                        continue
                    
                    parts = line.split()
                    # 确保数据列足够 (列0为基因1, 列2为基因2)
                    if len(parts) < 3:
                        continue

                    gene1 = parts[0]
                    gene2 = parts[2]

                    # 解析染色体编号
                    chr1 = parse_chromosome(gene1)
                    chr2 = parse_chromosome(gene2)

                    # 检查是否在目标集合中 (利用Set的O(1)查找)
                    if (chr1, chr2) in target_pairs:
                        # 写入结果 (制表符分隔)
                        f_out.write(f"{gene1}\t{gene2}\n")
                        filtered_count += 1
                    
                    line_processed += 1

            print(f"处理完成: {out_file}")
            print(f"统计: 扫描基因行 {line_processed} 条，提取符合条件 {filtered_count} 对")

        except Exception as e:
            print(f"[错误] 处理文件 {in_file} 时发生错误: {e}")

# -----------------------------------------------------------------------------
# 主程序入口
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="提取共线性文件中的特定染色体直系同源基因对",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # 命令行参数定义
    parser.add_argument(
        '-i', '--input', 
        required=True, 
        help='输入文件的匹配模式 (必须用引号包裹，防止Shell提前展开)\n例如: -i "*.block.rr.txt"'
    )
    parser.add_argument(
        '-o', '--output-suffix', 
        default='.pseu.ortologs', 
        help='输出文件的后缀名 (默认: .pseu.ortologs)'
    )
    parser.add_argument(
        '-c', '--chroms', 
        default='1,1;11,11', 
        help='需要筛选的染色体对，用分号分隔组，逗号分隔对 (默认: 1,1;11,11)\n例如: -c "1,1;2,2;3,5"'
    )

    args = parser.parse_args()

    # 1. 解析目标染色体
    target_pairs = parse_target_chroms(args.chroms)
    
    # 2. 执行处理
    process_block_files(args.input, args.output_suffix, target_pairs)
    
    print(f"--------------------------------------------------")
    print("所有任务已完成。")

if __name__ == "__main__":
    main()