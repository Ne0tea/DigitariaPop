'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-06-16 16:25:07
LastEditors: Ne0tea
LastEditTime: 2025-06-17 23:31:54
'''
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import zscore
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
def extract_chromosome(gene_id):
    match1 = re.search(r'Dexi([1-9][A-Za-z]?)', gene_id)
    match2 = re.search(r'(Chr\d+_RagTag)', gene_id)
    match3 = re.search(r'(Chr\d+)', gene_id)

    if match1:
        return match1.group(1)
    elif match2:
        return match2.group(1)
    elif match3:
        return match3.group(1)
    else:
        return None

def transform_row(row):
    row_name = row.name
    value_dic = {}
    for sp in row.index:
        cur_value = row[sp] / sp_ctg_gene_num[sp][row_name]
        value_dic[sp] = cur_value
    out = pd.Series(value_dic)
    return out


chr_trans_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\0_Digitaria_chr_trans'
ctg_homochr = {}
with open(chr_trans_file, 'r') as chr_trans_f:
    for i  in chr_trans_f:
        chr_name, subg_name, subg = i.split()[0], i.split()[1], i.split()[2]
        ctg_homochr[chr_name] = [subg_name, subg]

chr_num_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Digitaria_chr_pfam_number'
sp_ctg_gene_num = {}
with open(chr_num_file, 'r') as chr_num_f:
    for i  in chr_num_f:
        chr_name, gene_number = i.split()[1], int(i.split()[0])
        subg = ctg_homochr[chr_name][1]
        if 'RagTag' in chr_name and subg == 'SubC':
            sp = 'Drad'
        elif 'RagTag' in chr_name and subg == 'SubD':
            sp = 'Dmil-D'
        elif 'RagTag' in chr_name and subg == 'SubE':
            sp = 'Dmil-E'
        elif 'Chr' in chr_name and subg == 'SubC':
            sp = 'Dsan-C'
        elif 'Chr' in chr_name and subg == 'SubD':
            sp = 'Dsan-D'
        elif 'Chr' in chr_name and subg == 'SubE':
            sp = 'Dsan-E'
        elif subg == 'SubA':
            sp = 'Dexi-A'
        elif subg == 'SubB':
            sp = 'Dexi-B'
        if sp in sp_ctg_gene_num:
            sp_ctg_gene_num[sp][ctg_homochr[chr_name][0]] = gene_number
        else:
            sp_ctg_gene_num[sp] = {ctg_homochr[chr_name][0]: gene_number}
print(sp_ctg_gene_num)

target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\UDP_gene'
phy_pfam_df = pd.DataFrame()
all_pfam_chr_number_dic = {}
for file_name in os.listdir(target_directory):
    cur_pfam_file = os.path.join(target_directory,file_name)
    cur_sp_name = file_name.split('_')[0]
    with open(cur_pfam_file, 'r') as f:
        all_pfam_chr_number_dic[cur_sp_name] = {}
        for line in f:
            gene_name = line.strip()
            cur_chr = extract_chromosome(gene_name)
            homo_chr = ctg_homochr[cur_chr][0]
            if homo_chr in all_pfam_chr_number_dic[cur_sp_name]:
                all_pfam_chr_number_dic[cur_sp_name][homo_chr] += 1
            else:
                all_pfam_chr_number_dic[cur_sp_name][homo_chr] = 1

pfam_df = pd.DataFrame(all_pfam_chr_number_dic)

fisher_results = pd.DataFrame()
for column in pfam_df.columns:
    a = pfam_df.loc["Chr03", column]
    b = sp_ctg_gene_num[column]['Chr03'] - a
    # 计算该列的非 Chr03 基因数
    c = pfam_df[column].sum() - a
    # 计算既不在该列也不在 Chr03 的基因数
    d = (sum(sp_ctg_gene_num[column].values()) - pfam_df[column].sum()) - c

    contingency_table = [
        [a, b],
        [c, d]
    ]
    odds_ratio, p_value = stats.fisher_exact(contingency_table)
    fisher_results = fisher_results._append(pd.Series({
        "p_value": p_value,
        "odds_ratio": odds_ratio,
        "contingency_table": contingency_table
    },name=column))
enrich_df = pd.DataFrame(fisher_results)
_, p_corrected, _, _ = multipletests(enrich_df['p_value'], 
                                    alpha=0.05, 
                                    method='bonferroni')
enrich_df['Qvalue'] = p_corrected
pfam_df = pfam_df.apply(transform_row, axis=1)
pfam_df = pfam_df.apply(zscore, nan_policy='omit')

pfam_df = pfam_df.fillna(pfam_df.min().min())
pfam_df = pfam_df.loc[['Chr01','Chr02','Chr03', 'Chr04','Chr05','Chr06', 'Chr07','Chr08','Chr09']]
pfam_df = pfam_df.rename(columns={'Dsan-E':'EH', 'Dsan-D':'DH', 'Dsan-C':'CH',
                                  'Dmil-E':'ET', 'Dmil-D':'DT', 'Drad':'CD',
                                  'Dexi-A':'AT', 'Dexi-B':'BT'})
pfam_df = pfam_df[['AT', 'BT', 'CD', 'DT', 'ET', 'CH', 'DH', 'EH']]
pfam_df = pfam_df.T
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(pfam_df, annot=True, cmap='coolwarm', 
            linewidths=1, linecolor='white', ax=ax)

for i, row_name in enumerate(enrich_df.index):
    ax.annotate(
        text=f"{enrich_df.loc[row_name, 'p_value']:.2e}",
        xy=(1.2, i+0.5),
        xycoords=ax.get_yaxis_transform(),
        ha='left', va='center',
        fontsize=12
    )
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
plt.subplots_adjust(right=0.8)
# plt.show()
plt.savefig(r"E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\UDP_gene\UDPGTs.pdf")