'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-07-28 11:47:59
LastEditors: Ne0tea
LastEditTime: 2025-08-01 11:28:53
'''
from venn import venn
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

def main(DE_gene_file):
    DE_expr_df = pd.read_csv(DE_gene_file, header=0)
    DE_expr_df[['Time']] = DE_expr_df["Proj"].str.rsplit("_", n=1, expand=True)[[1]]
    DE_expr_df = DE_expr_df[(DE_expr_df['padj'] < 0.05) & (DE_expr_df['log2FoldChange'].apply(abs) >= 1)]
    DE_expr_df['Regulation'] = DE_expr_df['log2FoldChange'].apply(lambda x: 'Up' if x > 0 else 'Down')

    Down_expr_df = DE_expr_df[DE_expr_df['Regulation'] == 'Down'].groupby('Time')['gene_id'].nunique().reset_index()
    Up_expr_df = DE_expr_df[DE_expr_df['Regulation'] == 'Up'].groupby('Time')['gene_id'].nunique().reset_index()
    print(Up_expr_df)

    grouped = DE_expr_df[DE_expr_df['Regulation'] == 'Up'].groupby('Time')['gene_id'].apply(set).reset_index()
    up_plot_data = {}
    for idx, row in grouped.iterrows():
        cur_id = row['Time']
        up_plot_data[cur_id] = row['gene_id']
    grouped = DE_expr_df[DE_expr_df['Regulation'] == 'Down'].groupby('Time')['gene_id'].apply(set).reset_index()
    Down_plot_data = {}
    for idx, row in grouped.iterrows():
        cur_id = row['Time']
        Down_plot_data[cur_id] = row['gene_id']

    up_group_dic = {k: up_plot_data[k] for k in up_plot_data}
    down_group_dic = {k: Down_plot_data[k] for k in Down_plot_data}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))

    venn(up_group_dic, ax=ax1)
    venn(down_group_dic, ax=ax2)
    # sns.barplot(x='Var', y='gene_id', hue='Regulation', data=Regulation_expr_df, ax=ax3)
    # ax3.spines['top'].set_visible(False)
    # ax3.spines['right'].set_visible(False)

    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\RvS_DEG_stat.pdf')

if __name__ == "__main__":
    DE_gene_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat_RvS_RvS_DEGs.csv'

    main(DE_gene_file)