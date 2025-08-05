'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-19 15:43:00
LastEditors: Ne0tea
LastEditTime: 2025-03-19 16:02:37
'''
import pandas as pd
import sys
import re
pfam_df = pd.read_table(sys.argv[1], usecols=list(range(0, 10)), header=None)

target_gene_list = []
with open(sys.argv[2], 'r') as f:
    for i in f:
        line=i.strip().split()
        target_gene_list.append(line[0])

pfam_df.iloc[:, 0] = pfam_df.iloc[:, 0].apply(lambda x: re.sub(r'Dsan-[CDE].', '', re.sub(r'.mRNA1', '', x.replace('_','.'))))
target_pfam = pfam_df[pfam_df.iloc[:, 0].isin(target_gene_list)]
target_pfam.to_csv(sys.argv[3], sep='\t', index=False, header=False)
