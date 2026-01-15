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

## 快速使用范例

### 步骤 1：筛选直系同源基因

**目标**：从共线性 Block 文件中提取特定染色体（如1号和11号）的直系同源对。

Bash

```
python 1.filterOrthlogs.py \
  -i "Lso_Lma.block.rr.txt" \
  -c "1,1;11,11"
```

> **结果**：生成文件 `Lso_Lma.block.rr.txt.pseu.ortologs`

------

### 步骤 2：构建四联子列表

**目标**：结合上一步生成的直系同源文件和旁系同源文件，组装四联子。

Bash

```
python 2.extractGeneQuartets.py \
  -o "Lso_Lma.block.rr.txt.pseu.ortologs" \
  -p "Lso.v.Lso.paralog" \
  -out "Lso_Lma.quartet"
```

> **结果**：生成文件 `Lso_Lma.quartet`

------

### 步骤 3：检测基因置换

**目标**：计算 Ka/Ks 并输出最终判定结果。 *(注：这里假设您的 CDS 序列文件名为 `Lso.cds.fasta` 和 `Lma.cds.fasta`，请根据实际文件名修改 -a 和 -b 参数)*

Bash

```
python 3.detetConver.py \
  -q "Lso_Lma.quartet" \
  -a "Lso.cds.fasta" \
  -b "Lma.cds.fasta" \
  -o "Lso_Lma.P.CV.PaPs.txt" \
  --boot 100
```

> **结果**：生成文件 `Lso_Lma.P.CV.PaPs.txt`



## 🚀 详细使用指南与示例

### 第一步：提取直系同源基因 (Filter Orthologs)

脚本: 1.filterOrthlogs.py

功能: 从共线性分析结果（如 ColinearrScan或MCScanX 输出的 .block 文件）中，提取特定染色体对的直系同源基因。

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


------

## 📊 结果解读 (Interpretation)

### 1. 📄 结果文件每列含义详解

该文件每一行代表一个“基因四联子”（Quartet）的分析结果。

| **列名**          | **含义说明**                                                 | **关键作用**                                                 |
| ----------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **QuartetID**     | 四联子编号，格式通常为 `物种A旁系1-物种A旁系2`。             | 唯一标识一组待检测的基因。                                   |
| **Ka_P1**         | 物种1（如 Lso）内两个**旁系基因**间的非同义突变率。          | 辅助参考，用于看蛋白差异。                                   |
| **Ks_P1**         | **物种1内两个旁系基因间的同义突变率**。                      | **核心指标**。代表旁系基因的分化时间。数值越小，说明越亲近。 |
| **Ka_P2 / Ks_P2** | 物种2（如 Lma）内旁系基因的 Ka/Ks 值。                       | 用于检测物种2是否发生置换。                                  |
| **Ka_O1**         | 第一对**直系同源基因**（如 Lso1 vs Lma1）的非同义突变率。    | 辅助参考。                                                   |
| **Ks_O1**         | **第一对直系同源基因间的同义突变率**。                       | **基准线（时间尺）**。代表两个物种的分化时间。               |
| **Ka_O2 / Ks_O2** | 第二对直系同源基因（如 Lso2 vs Lma2）的 Ka/Ks 值。           | 辅助基准线，用于验证 Ks_O1 的可靠性。                        |
| **Conv_Sp1**      | **物种1是否发生置换的判定结果** (`Y`/`N`)。 判断逻辑：若 `Ks_P1 < Ks_O1`，则为 `Y`。 | 直接告诉你是否检测到异常。                                   |
| **Conv_Sp2**      | 物种2是否发生置换的判定结果 (`Y`/`N`)。                      | 同上。                                                       |
| **Boot_Prob_Sp1** | **物种1置换结果的 Bootstrap 支持率** (0.00~1.00)。           | **可信度指标**。数值越高（如 >0.9），结果越可靠。            |
| **Boot_Prob_Sp2** | 物种2置换结果的 Bootstrap 支持率。                           | 同上。                                                       |

------

### 2. 🔍 如何区分完全基因置换 (WCV) 与部分基因置换 (PCV)

虽然脚本输出的是整条序列的平均值，但结合 **Ks 数值** 和 **Bootstrap 支持率** 的特征，我们可以对置换类型进行分类推断。

#### ✅ 完全基因置换 (Whole Gene Conversion, WCV)

**定义**：基因 A 将基因 B 的**整个序列**完全覆盖，导致两者全长几乎一模一样。

- **判定标准（需同时满足）**：
  1. **判定列 (Conv)**：显示为 **`Y`**。
  2. **Ks 值 (Ks_P1)**：**极低，接近于 0**（例如 `< 0.01` 或显示为 `0.0000`）。这说明序列差异极小。
  3. **支持率 (Boot_Prob)**：**极高**（通常 **`> 0.90`** 甚至 `1.00`）。
     - *原理*：因为整条序列都发生了置换，信号非常强且一致，无论 Bootstrap 怎么随机抽样，都能检测到置换信号。

#### ⚠️ 部分基因置换 (Partial Gene Conversion, PCV) 候选

**定义**：基因 A 只覆盖了基因 B 的**一段区域**（如某个外显子或功能域），形成“镶嵌”结构。

- **判定标准（特征组合）**：
  1. **判定列 (Conv)**：显示为 **`Y`**。
  2. **Ks 值 (Ks_P1)**：**显著低于基准线 (`Ks_O1`)，但明显大于 0**（例如 `Ks_O1=0.4`，而 `Ks_P1=0.15`）。
     - *原理*：置换区域的 Ks 接近 0，但未置换区域的 Ks 很高（正常进化），两者一平均，就得到了一个中间值。
  3. **支持率 (Boot_Prob)**：**中等水平**（通常在 **`0.50 - 0.90`** 之间）。
     - *原理*：Bootstrap 重采样时，如果抽到了置换区域的位点，结果就是 Y；如果抽到了未置换区域，结果就是 N。信号不稳定导致支持率下降。
     - **建议**: 对于这类基因，建议后续进行 **滑动窗口分析 (Sliding Window Analysis)** 以确定具体的置换断点。

### 总结：操作指南

您可以根据以下逻辑对 `Lso_Lma.P.CV.PaPs.txt` 进行分类：

| **类别**                             | **筛选条件 (基于文章 WCV-I 方法)**          | **生物学含义**                                               |
| ------------------------------------ | ------------------------------------------- | ------------------------------------------------------------ |
| **完全基因置换 (Whole)**             | `Conv = Y` **且** `Boot_Prob ≥ 0.90`        | 极高置信度。整条基因序列都极度相似，甚至 Ks ≈ 0。            |
| **潜在部分置换 (Partial Candidate)** | `Conv = Y` **且** `0.50 ≤ Boot_Prob < 0.90` | 信号存在但不稳定。暗示置换可能只发生在基因的某个区域，导致全长支持率不高。 |
| **无置换 / 噪音**                    | `Conv = N` **或** `Boot_Prob < 0.50`        | 正常进化，或者信号太弱无法区分是噪音还是微小置换。           |



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



