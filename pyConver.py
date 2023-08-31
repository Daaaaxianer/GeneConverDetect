# -*- coding:utf-8 -*-
# @FileName :pyConver.py
# @Time     :2023/4/18 22:42
# @Author   :Xian

import pandas as pd

# 设置pandas显示所有列和行
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
# 设置pandas字符宽度为100
pd.set_option('display.width', 200)

### 1. 读取共线性alignment，提取同源基因四联子(quartet)
orthal = pd.read_csv("./alignment/Osa.Ssp.alignment.csv", header=None, names=["os", "ss1", "ss2", "ss3", "ss4"], na_values=["."], usecols=[0, 1, 2, 3, 4])
refal = pd.read_csv('./alignment/Osa.Osa.alignment.csv', header=None, names=['os1', 'os2'], na_values=['.'])
nrefal = pd.read_csv('./alignment/Ssp.Ssp.alignment.csv', header=None, names=['ss1', 'ss2', 'ss3', 'ss4'], na_values=['.'], usecols=[0, 1, 2, 3])

print(orthal.head(10))
print(refal.head(250))
# print(nrefal.head(10))
print("#" * 100)

# 利用orthal的'os'为索引，为refal增加三列
al_df = pd.merge(refal, orthal, left_on='os1', right_on='os')
# print(al_df.head(250))
print("#" * 100)
al_df = pd.merge(al_df, orthal, left_on='os2', right_on='os', how='left', suffixes=['_1', '_2']).drop(['os1', 'os2'], axis=1)
print(al_df.head(250))

# 提取四联子
# al_df = al_df.loc[:, ['ss1', 'ss2', 'ss3', 'ss4']]

### 2. 对四联子进行多序列比对，删除gap > 50% 或 identity < 40% 的四联子


### 3. 对剩余四联子构建基因树


### 4. 比较置换基因对的ks值

