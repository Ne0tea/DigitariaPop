'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-06-20 21:14:11
LastEditors: Ne0tea
LastEditTime: 2025-08-02 14:37:33
'''
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from joypy import joyplot
import numpy as np
from scipy import stats
import seaborn as sns
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
ignore_material = ['YZGJ2_1','NT3_1']
joy_color = {'Resistant':'#dda15e','Sensitive':'#ccd5ae'}
Doamin_color = {'ALS':'#ffddd2', 'EPSPS':'#edf6f9', 'UGT':'#e29578', 'AKR':'#f9ebc7', 'ABC':'#ccd5ae', 'NB-ARC':'#d4a373', 'GST':'#83c5be', 'P450':'#006d77'}
sub_color_set={'SubC':'#003049','SubD':'#d62828','SubE':'#f77f00'}

def is_normal_with_mean_1(data, alpha=0.05):
    # 1. 正态性检验(Shapiro-Wilk检验)
    _, p_value_normal = stats.shapiro(data)
    is_normal = p_value_normal > alpha
    t_stat, p_value_mean = stats.normaltest(data)
    print(np.mean(data))
    is_mean_1 = p_value_mean > alpha

    return p_value_normal, p_value_mean

def draw_boxplot(all_gene_depth_df):
    # plt.figure(figsize=(3, 5))
    groups = all_gene_depth_df['Domain'].unique()
    all_gene_depth_df = all_gene_depth_df[all_gene_depth_df['CN'] < 4]
    fig, axes = plt.subplots(8, 1, figsize=(4, 5), sharex=True)

    for i, group in enumerate(groups):
        ax = axes[i]
        sns.histplot(
            data=all_gene_depth_df[all_gene_depth_df['Domain'] == group],
            x='CN',
            ax=ax,
            kde=True,
            bins=100,linewidth=0,
            # color=None,
            color=Doamin_color[group],
            alpha=0.7,
            label=f'{group}'
        )
        ax.set_xlim(0,4)
        ax.set_ylabel(group, rotation=0, ha='right', va='center')
        ax.axvline(x=1,color='darkred',linestyle='--',linewidth=2)
        ax.legend().set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_yticks([])

    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\Domain_CN.pdf', format='pdf')
    # sns.boxplot(x="Domain", y="CN", data=all_gene_depth_df, palette=Doamin_color, showfliers=False)
    # sns.stripplot(x="Domain", y="CN", data=all_gene_depth_df, palette=Doamin_color, alpha=0.3)

    # plt.savefig(os.path.join(out_dir,prefix+'_'+str(i*100)+'.pdf'), format='pdf')
    # plt.close()

def main(In_article_sample, ratio_file_dir, herbicide_resistance_file):
    article_sample_list = []
    with open(In_article_sample, 'r') as article_samplef:
        for i in article_samplef:
            line = i.strip().split()
            article_sample_list.append(line[0])

    HB_RR_Df = pd.read_table(herbicide_resistance_file, sep="\t").set_index("Vcf_id")
    HB_RR_Df = HB_RR_Df[HB_RR_Df['Cluster_in_material'] == 'C5']
    C5_material_list = HB_RR_Df[(HB_RR_Df['Cluster_in_material'] == 'C5') & (HB_RR_Df['Type'] == 'Resistant')].index.to_list()

    all_gene_depth_df = pd.DataFrame()
    for i in os.listdir(ratio_file_dir):
        # if i.endswith('_ratio.txt') and 'uniq' not in i:
        if i.endswith('uniq_ratio.txt'):
            cur_file_dir = os.path.join(ratio_file_dir,i)
            prefix = i.split('_')[0]
            cur_depth_df = pd.read_table(cur_file_dir)
            cur_depth_df.columns = ["ID"] + cur_depth_df.columns[1:].tolist()
            cur_depth_df = cur_depth_df.set_index('ID')
            cur_depth_df = cur_depth_df.T.reset_index(names="ID")

            # if 'ALS' in i:
            #     # article_sample_list.remove('NT3_1')
            #     # article_sample_list.remove('YZGJ2_1')
            #     p_normal, p_mean = is_normal_with_mean_1(cur_depth_df.loc[cur_depth_df['ID'].isin(article_sample_list), 'Chr16.1966'])
            #     print(p_normal, p_mean)
            #     p_normal, p_mean = is_normal_with_mean_1(cur_depth_df.loc[cur_depth_df['ID'].isin(article_sample_list), 'Chr17.1896'])
            #     print(p_normal, p_mean)
            #     p_normal, p_mean = is_normal_with_mean_1(cur_depth_df.loc[cur_depth_df['ID'].isin(article_sample_list), 'Chr18.2017'])
            #     print(p_normal, p_mean)
            #     continue
            # else:
            #     continue

            df_long = pd.melt(cur_depth_df,
                                id_vars=['ID'],
                                value_vars=cur_depth_df.columns.tolist(),
                                var_name='GeneID',
                                value_name='CN'
                            )
            df_long['Domain'] = prefix
            all_gene_depth_df = pd.concat([all_gene_depth_df, df_long])

    # all_gene_depth_df = all_gene_depth_df[all_gene_depth_df['ID'].isin(C5_material_list)]
    all_gene_depth_df = all_gene_depth_df[~all_gene_depth_df['ID'].isin(ignore_material)]
    all_gene_depth_df = all_gene_depth_df[all_gene_depth_df['ID'].isin(article_sample_list + ['Type'])]
    all_gene_depth_df = all_gene_depth_df[all_gene_depth_df['CN'] > 3]
    print(all_gene_depth_df.groupby(['Domain'])['ID'].nunique())
    # all_gene_depth_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\NTSR_CNV_mean.csv')
    # xx = all_gene_depth_df.groupby(['Domain', 'ID'])['CN'].mean().sort_values(ascending=False)
    # xx.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\NTSR_CNV_mean.csv')
    # draw_boxplot(all_gene_depth_df)

if __name__ == "__main__":
    In_article_sample = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\Article_sample.list'
    ratio_file_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis'
    herbicide_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'

    main(In_article_sample, ratio_file_dir, herbicide_resistance_file)