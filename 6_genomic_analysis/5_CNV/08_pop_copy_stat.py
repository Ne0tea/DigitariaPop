'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-18 23:26:35
LastEditors: Ne0tea
LastEditTime: 2025-07-27 23:30:15
'''
'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-15 23:59:34
LastEditors: Ne0tea
LastEditTime: 2025-06-20 20:55:20
'''
import sys
import os
import pandas as pd
import PyComplexHeatmap as pch
import numpy as np
from joypy import joyplot
import matplotlib.pyplot as plt
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
ignore_material = ['YZGJ2_1','NT3_1']
joy_color = {'Resistant':'#dda15e','Sensitive':'#ccd5ae'}
sub_color_set={'SubC':'#003049','SubD':'#d62828','SubE':'#f77f00'}

def draw_joyplot(depth_df, Vcf_RR_dic):
    depth_plot_df = pd.melt(
    depth_df,
    id_vars=['geneID'],
    var_name='Vcf_id',  # 新列名（存储转换前的列名）
    value_name='Depth'  # 新列名（存储转换前的值）
    )
    depth_plot_df['Type'] = depth_plot_df['Vcf_id'].map(Vcf_RR_dic)
    depth_plot_df.dropna(inplace=True)
    print(depth_plot_df)
    depth_plot_df = depth_plot_df.pivot(
        index=['geneID', 'Vcf_id'],  # 保留不变的列
        columns='Type',     # 转换为新列的列
        values='Depth'         # 新列的值
    ).reset_index()

    fig,ax = joyplot(depth_plot_df, by='geneID', column=['Resistant', 'Sensitive'], figsize= (4,8), legend=True,hist=True,bins=20,
            color=['#dda15e','#ccd5ae'], linewidth=.8, grid=True, alpha=.3, x_range=[-1, 5])

    plt.subplots_adjust(top=0.95,bottom=0.1)
    plt.show()
    # plt.savefig(os.path.join(out_dir,prefix+'_'+str(i*100)+'.pdf'), format='pdf')
    plt.close()

def main(In_article_sample, ratio_file_dir, herbicide_resistance_file):
    article_sample_list = []
    with open(In_article_sample, 'r') as article_samplef:
        for i in article_samplef:
            line = i.strip().split()
            article_sample_list.append(line[0])

    all_gene_depth_df = pd.DataFrame()
    for i in os.listdir(ratio_file_dir):
        if i.endswith('_ratio.txt'):
            cur_file_dir = os.path.join(ratio_file_dir,i)
            prefix = i.split('_')[0]
            cur_depth_df = pd.read_table(cur_file_dir)
            cur_depth_df.columns = ["geneID"] + cur_depth_df.columns[1:].tolist()
            cur_depth_df = cur_depth_df.set_index('geneID')
            cur_depth_df = cur_depth_df.T
            cur_depth_df.loc['Type'] = prefix
            # print(cur_depth_df)
            # cur_depth_df['Max_'+prefix] = cur_depth_df.mean(axis=1)
            # all_gene_depth_df = pd.concat([all_gene_depth_df,cur_depth_df['Max_'+prefix]],axis=1)
            all_gene_depth_df = pd.concat([all_gene_depth_df,cur_depth_df],axis=1)
    all_gene_depth_df = all_gene_depth_df[~all_gene_depth_df.index.isin(ignore_material)]
    all_gene_depth_df = all_gene_depth_df[all_gene_depth_df.index.isin(article_sample_list + ['Type'])]
    all_gene_depth_df = all_gene_depth_df.loc[:, (all_gene_depth_df.iloc[:-1] > 3).any()]
    print(all_gene_depth_df)
    all_gene_depth_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\NTSR_CNV.csv')

    # HB_RR_Df = pd.read_table(herbicide_resistance_file, sep="\t").set_index("Vcf_id")
    # HB_RR_Df = HB_RR_Df[HB_RR_Df['Cluster_in_material'] == 'C5']
    # C5_material_list = HB_RR_Df[HB_RR_Df['Cluster_in_material'] == 'C5'].index.to_list()

    # Vcf_RR_dic = HB_RR_Df['Type'].to_dict()
    # all_gene_depth_df = all_gene_depth_df[all_gene_depth_df.index.isin(C5_material_list)]
    # # all_gene_depth_df = np.log2(all_gene_depth_df)
    # type_expr_df = pd.DataFrame(Vcf_RR_dic, index=[0]).T
    # type_expr_df.columns = ['Type']
    # # print(type_expr_df)

if __name__ == "__main__":
    In_article_sample = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\Article_sample.list'
    ratio_file_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis'
    herbicide_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
    # target_gene_file = []
    main(In_article_sample, ratio_file_dir, herbicide_resistance_file)