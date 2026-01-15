# Gene Conversion Detection Pipeline

# 基因置换检测分析流程

本流程包含三个 Python 脚本，旨在从全基因组复制（WGD）或多倍体背景下，系统性地筛选并检测 **基因置换（Gene Conversion）** 事件。

流程基于 **四联子（Quartet）** 模型，结合 **共线性分析**、**同义突变率（Ks）计算** 以及 **Bootstrap 统计检验**，能够识别出违反正常进化时钟的异常基因对。

------

## 🛠️ 环境与依赖 (Requirements)

在运行之前，请确保您的环境满足以下要求：

### 1. 软件依赖

本流程的核心计算依赖 **ClustalW** 进行多序列比对，请务必安装并添加到环境变量 PATH 中。

- **Linux (Ubuntu/Debian)**: `sudo apt-get install clustalw`
- **macOS**: `brew install clustal-w`
- **Conda**: `conda install -c bioconda clustalw`
- **验证安装**: 在终端输入 `clustalw2 -help`，应能看到帮助信息。

### 2. Python 库

Bash

```
pip install biopython argparse
```

------

## 📂 脚本列表与功能

| **顺序** | **脚本名称**               | **功能描述**                                                 |
| -------- | -------------------------- | ------------------------------------------------------------ |
| **01**   | `1.filterOrthlogs.py`      | **筛选直系同源**：从共线性 Block 文件中提取特定染色体对的基因对。 |
| **02**   | `2.extractGeneQuartets.py` | **构建四联子**：结合旁系同源列表，组装 `(Para1, Para2)-(Ortho1, Ortho2)` 结构。 |
| **03**   | `3.detetConver.py`         | **检测置换**：进行序列比对、计算 Ka/Ks、判定置换并执行 Bootstrap 验证。 |

------

## 🚀 详细使用指南与示例

### 第一步：提取直系同源基因 (Filter Orthologs)

脚本: 1.filterOrthlogs.py

功能: 从共线性分析结果（如 MCScanX 输出的 .block 文件）中，提取特定染色体对的直系同源基因。

#### 📖 参数说明

- `-i`, `--input`: 输入文件路径（必填，支持通配符 `*`，如 `"*.txt"`，需加引号）。
- `-o`, `--output-suffix`: 输出文件的后缀（默认 `.pseu.ortologs`）。
- `-c`, `--chroms`: 目标染色体对。格式为 `chrA,chrB;chrC,chrD`。默认 `1,1;11,11`。

#### 💡 使用示例

**示例 1：基础用法（默认筛选 Chr1-Chr1 和 Chr11-Chr11）**

Bash

```
python 1.filterOrthlogs.py -i "data/Lso_Lma.block.rr.txt"
```

> *结果生成：`data/Lso_Lma.block.rr.txt.pseu.ortologs`*

**示例 2：自定义筛选特定的染色体对（如 2号对2号，3号对5号）**

Bash

```
python 1.filterOrthlogs.py \
  -i "data/Lso_Lma.block.rr.txt" \
  -c "2,2;3,5"
```

**示例 3：批量处理当前目录下所有的 block 文件**

Bash

```
python 1.filterOrthlogs.py \
  -i "*.block.rr.txt" \
  -o ".filtered.ortho.txt"
```

------

### 第二步：构建四联子列表 (Extract Quartets)

脚本: 2.extractGeneQuartets.py

功能: 利用上一步生成的直系同源表，结合物种内部的旁系同源表，构建用于检测的四基因组合。

#### 📖 参数说明

- `-o`, `--ortho`: 第一步生成的直系同源文件路径（必填）。
- `-p`, `--para`: 旁系同源基因列表文件（必填，格式至少包含两列基因ID）。
- `-out`, `--output`: 输出的四联子列表文件路径（默认 `Lso_Lma.quartet`）。

#### 💡 使用示例

**示例 1：使用默认文件名快速运行**

Bash

```
# 前提：当前目录下存在 Lso_Lma.block.rr.txt.pseu.ortologs 和 Lso.v.Lso.paralog
python 2.extractGeneQuartets.py
```

**示例 2：指定自定义的输入输出文件**

Bash

```
python 2.extractGeneQuartets.py \
  -o "data/my_orthologs.txt" \
  -p "data/my_paralogs.txt" \
  -out "results/my_quartet_list.txt"
```

**示例 3：处理不同物种的数据**

Bash

```
python 2.extractGeneQuartets.py \
  -o "Rice_Maize.orthologs" \
  -p "Rice.paralogs" \
  -out "Rice_Maize.quartet"
```

------

### 第三步：检测基因置换 (Detect Conversion)

脚本: 3.detetConver.py

功能: 这是核心分析步骤。脚本会自动提取序列、翻译蛋白、调用 ClustalW 比对、回译 DNA 并计算 Ka/Ks。

#### 📖 参数说明

- `-q`, `--quartet`: 第二步生成的四联子文件（必填）。
- `-a`, `--fasta1`: 物种 1 的 CDS 序列文件 (FASTA格式，必填)。
- `-b`, `--fasta2`: 物种 2 的 CDS 序列文件 (FASTA格式，必填)。
- `-o`, `--output`: 最终结果输出文件路径（必填）。
- `--boot`: Bootstrap 重采样次数（默认 100）。建议设为 1000 以获得出版级可信度。

#### 💡 使用示例

**示例 1：快速测试（默认 Bootstrap 100 次）**

Bash

```
python 3.detetConver.py \
  -q "results/Lso_Lma.quartet" \
  -a "genomes/Lso.cds.fasta" \
  -b "genomes/Lma.cds.fasta" \
  -o "results/final_report.txt"
```

**示例 2：高精度分析（提高 Bootstrap 次数至 1000）**

Bash

```
python 3.detetConver.py \
  -q "results/Lso_Lma.quartet" \
  -a "genomes/Lso.cds.fasta" \
  -b "genomes/Lma.cds.fasta" \
  -o "results/publication_quality_report.txt" \
  --boot 1000
```

**示例 3：仅计算 Ka/Ks 值（设 Bootstrap 为 0，速度最快，用于初筛）**

Bash

```
python 3.detetConver.py \
  -q "results/large_dataset.quartet" \
  -a "genomes/Lso.cds.fasta" \
  -b "genomes/Lma.cds.fasta" \
  -o "results/screening_report.txt" \
  --boot 0
```

------

## 📊 结果解读 (Interpretation)

结果文件 (`final_report.txt`) 包含多列数据，以下是关键列的生物学含义及**完全/部分置换**的区分方法。

### 1. 关键列说明

| **列名**          | **含义**                                                     |
| ----------------- | ------------------------------------------------------------ |
| **QuartetID**     | 四联子编号 (Para1-Para2)                                     |
| **Ks_P1**         | 物种1内 **旁系基因** 间的 Ks 值（同义突变率）。              |
| **Ks_O1**         | 物种间 **直系基因** 间的 Ks 值（基准时间尺）。               |
| **Conv_Sp1**      | **判定结果**。若 `Ks_P1 < Ks_O1`，则标记为 `Y` (Yes)，否则为 `N` (No)。 |
| **Boot_Prob_Sp1** | **Bootstrap 支持率** (0.00 - 1.00)。表示随机重采样验证中，支持置换结论的比例。 |

### 2. 判读逻辑：正常 vs 置换

- **正常进化 (N)**:
  - `Ks_P1` (大) > `Ks_O1` (小)。
  - 说明旁系基因分家很早，积累了大量差异。
- **基因置换 (Y)**:
  - `Ks_P1` (小) < `Ks_O1` (大)。
  - 说明旁系基因之间发生了序列交换，导致它们“返老还童”，看起来比刚分家的物种还亲近。

### 3. 进阶判读：完全置换 vs 部分置换

虽然脚本基于全长序列给出一个判读结果，但结合 `Ks` 数值和 `Bootstrap` 可以区分置换类型：

#### A. 完全基因置换 (Whole Gene Conversion)

两个基因在全长范围内几乎一模一样。

- **特征**:
  - `Conv_Sp1` = **Y**
  - `Ks_P1` **极低** (接近 0，例如 < 0.005)。
  - `Boot_Prob` **极高** (通常 > 0.95 或 1.00)。
- **含义**: 基因 A 把 基因 B 彻底覆盖了。

#### B. 部分基因置换 (Partial Gene Conversion)

基因呈“马赛克”结构，只有一段发生了置换，其余部分差异很大。

- **特征**:
  - `Conv_Sp1` = **Y** (或 N，视片段大小而定)。
  - `Ks_P1` **偏低**，但不是 0 (处于中间值，例如 0.05)。这是因为置换区的低 Ks 被未置换区的高 Ks 稀释了。
  - `Boot_Prob` **中等** (例如 0.50 - 0.90)。因为随机抽样时，有时抽到置换区（支持Y），有时抽到未置换区（支持N），导致支持率震荡。
- **建议**: 对于这类基因，建议后续进行 **滑动窗口分析 (Sliding Window Analysis)** 以确定具体的置换断点。

------

## ⚠️ 常见问题 (Troubleshooting)

1. **报错 `Command 'clustalw2' not found`**:
   - 请检查 ClustalW 是否安装。如果已安装，请确保它在系统的 PATH 中，或者修改 `3.detetConver.py` 中的 `cmd` 变量指向绝对路径。
2. **结果中全是 `N`**:
   - 检查 `-a` 和 `-b` 的 FASTA 文件中 ID 是否与四联子文件中的 ID 完全匹配。
   - 基因置换是相对罕见的事件，如果没有是正常的。
3. **Bootstrap 运行很慢**:
   - 这是正常的。ClustalW 需要对每次重采样进行比对。如果数据量大，可以将 `--boot` 设为 0 先快速筛选候选者，再对候选者进行高次数 Bootstrap。

------

## 📜 引用与协议

本流程脚本遵循 MIT 协议。如果您在研究中使用了此流程，请保留脚本头部说明。