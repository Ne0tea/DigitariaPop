'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-27 15:41:33
LastEditors: Ne0tea
LastEditTime: 2025-07-27 16:32:12
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat 2.txt'
R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_count_matrix.csv'
S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_count_matrix.csv'

def read_stringtie_out(file):
    result_df = pd.read_csv(file, header=0)
    result_df = result_df.set_index('gene_id')
    return result_df

def make_treat_dic(treat_file):
    treat_name_SRR={}
    CK_dic = {}
    with open(treat_file, 'r') as gf:
        for i in gf:
            line=i.strip().split()
            srr_id = line[0]
            library_type = line[1].split('|')[1]
            treat_name = line[1].split('|')[0]
            # if srr_id == '15-2-0-3': continue
            if library_type == 'CK': 
                if treat_name.split('_')[0] in CK_dic:
                    CK_dic[treat_name.split('_')[0]].append(srr_id)
                else:
                    CK_dic[treat_name.split('_')[0]] = [srr_id]
                continue
            if treat_name in treat_name_SRR:
                treat_name_SRR[treat_name][library_type].append(srr_id)
            else:
                treat_name_SRR[treat_name]={'CK':[srr_id],'treatment':[]} if library_type == 'CK' else {'CK':[],'treatment':[srr_id]}
    for i in treat_name_SRR:
        cur_type = i.split('_')[0]
        treat_name_SRR[i]['CK'] = CK_dic[cur_type]

    return treat_name_SRR

def calculate_intragroup_correlation(df, treatment):
    """计算单个处理组内重复样本的基因表达相关性"""
    group_samples = [col for col in df.columns if col.startswith(treatment)]
    return df[group_samples].corr(method='spearman')

R_count_df = read_stringtie_out(R_count_file)
S_count_df = read_stringtie_out(S_count_file)
# S_count_df = S_count_df.drop(columns=['15-2-0-3'])
expr_matrix = pd.merge(R_count_df, S_count_df, how='outer', left_index=True, right_index=True)
# treat_name_SRR = make_treat_dic(treat_file)
treat_name = [ x.rsplit('-',1)[0] for x in list(expr_matrix.columns)]
treat_name = list(set(treat_name))
treat_name = [x for x in treat_name if 'Y' in x or '0' in x]

corr_results = {
    treatment: calculate_intragroup_correlation(expr_matrix, treatment)
    for treatment in treat_name  # 自动识别所有处理组
}

fig, axes = plt.subplots(4, 4, 
                        figsize=(8, 8),
                        sharey=False)

for idx, (treatment, corr) in enumerate(corr_results.items()):
    row = idx // 4
    col = idx % 4
    ax = axes[row, col]

    ax.set_axisbelow(True)
    ax.grid(True, color='lightgray', linestyle='--', linewidth=1, zorder=0)
    group_samples = [s for s in treat_name if s.startswith(treatment)]
    subset_corr = pd.DataFrame(corr)
    mask = np.triu(np.ones_like(subset_corr, dtype=bool))

    subset_corr.columns = list(map(lambda x: x.replace('15-2','S'), list(subset_corr.columns)))
    subset_corr.index = list(map(lambda x: x.replace('15-2','S'), list(subset_corr.index)))
    subset_corr.columns = list(map(lambda x: x.replace('21-17','R'), list(subset_corr.columns)))
    subset_corr.index = list(map(lambda x: x.replace('21-17','R'), list(subset_corr.index)))
    sns.heatmap(subset_corr, 
                ax=ax,
                cmap='coolwarm',
                vmin=-1, 
                vmax=1,
                mask=mask,
                annot=True, 
                fmt=".2f",
                cbar=False,
                square=True, zorder=1
                )

    ax.set_xticklabels(list(subset_corr.columns), rotation=45, ha='right')
    ax.set_yticklabels(list(subset_corr.index), rotation=0)

axes[3][2].set_facecolor('none')
axes[3][2].axis('off')
axes[3][3].set_facecolor('none')
axes[3][3].axis('off')
# 添加共享的颜色条
fig.colorbar(axes[0][0].collections[0], ax=axes[3,2],
             location='left',
            #  shrink=0.6,
             label='Spearman Correlation')

plt.tight_layout()
# plt.show()
plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Supplementary Fig. Pearson correlation of different library.pdf', format='pdf')