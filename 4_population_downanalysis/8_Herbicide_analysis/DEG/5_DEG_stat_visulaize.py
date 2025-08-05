'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-07-28 11:47:59
LastEditors: Ne0tea
LastEditTime: 2025-08-01 10:56:53
'''
from venn import venn
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

def make_treat_dic(treat_file):
    treat_name_SRR={}
    CK_list = []
    with open(treat_file, 'r') as gf:
        for i in gf:
            line=i.strip().split()
            srr_id = line[0]
            libary_type = line[1].split('|')[1]
            if libary_type == 'CK':
                CK_list.append(srr_id)
                continue
            treat_name = line[1].split('|')[0].split('_')[1]
            if srr_id == '15-2-0-3': continue
            if treat_name in treat_name_SRR:
                treat_name_SRR[treat_name].append(srr_id)
            else:
                treat_name_SRR[treat_name]=[srr_id]
    return treat_name_SRR, CK_list

def main(treat_file, DE_gene_file, R_count_file, S_count_file):
    treat_name_SRR, CK_list = make_treat_dic(treat_file)
    DE_expr_df = pd.read_table(DE_gene_file, header=0)
    DE_expr_df[['Var', 'Ht', 'Time']] = DE_expr_df["Proj"].str.split("_", expand=True)[[0,1,2]]
    DE_expr_df = DE_expr_df[(DE_expr_df['padj'] < 0.05) & (DE_expr_df['log2FoldChange'].apply(abs) >= 1)]
    DE_expr_df['Regulation'] = DE_expr_df['log2FoldChange'].apply(lambda x: 'Up' if x > 0 else 'Down')

    Regulation_expr_df = DE_expr_df.groupby(['Var', 'Regulation'])['gene_id'].nunique().reset_index()

    grouped = DE_expr_df.groupby(['Var', 'Time'])['gene_id'].apply(set).reset_index()
    plot_data = {}
    for idx, row in grouped.iterrows():
        cur_id = row['Var'] + '_' + row['Time']
        plot_data[cur_id] = row['gene_id']

    early_group_dic = {k: plot_data[k] for k in plot_data if k in ['Resistance_6h', 'Sensitive_6h']}
    late_group_dic = {k: plot_data[k] for k in plot_data if k in ['Resistance_24h', 'Sensitive_24h']}
    print(early_group_dic)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 3))

    venn(early_group_dic, ax=ax1)
    venn(late_group_dic, ax=ax2)
    sns.barplot(x='Var', y='gene_id', hue='Regulation', data=Regulation_expr_df, ax=ax3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    plt.show()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\DEG_stat.pdf')

if __name__ == "__main__":
    treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat.txt'
    DE_gene_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat_dif_gene.txt'
    R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_fpkm_matrix.csv'
    S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_fpkm_matrix.csv'

    main(treat_file, DE_gene_file, R_count_file, S_count_file)