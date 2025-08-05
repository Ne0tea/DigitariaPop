'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-15 23:59:34
LastEditors: Ne0tea
LastEditTime: 2025-08-02 14:11:32
'''
import pandas as pd
import seaborn as sns
from joypy import joyplot
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
ignore_material = ['YZGJ2_1','NT3_1']
joy_color = {'Resistant':'#dda15e','Sensitive':'#ccd5ae'}
sub_color_set={'SubC':'#003049','SubD':'#d62828','SubE':'#f77f00'}
Doamin_color = {'ALS':'#ffddd2', 'EPSPS':'#edf6f9', 'UGT':'#e29578', 'AKR':'#f9ebc7', 'ABC':'#ccd5ae', 'NB-ARC':'#d4a373', 'GST':'#83c5be', 'P450':'#006d77'}

def plot_custom_joint(ax_main, ax_y, depth_df, Vcf_RR_dic, domain):
    depth_plot_df = pd.melt(
    depth_df,
    id_vars=['geneID'],
    var_name='Vcf_id',
    value_name='Depth'
    )
    plot_df= pd.merge(depth_plot_df, Vcf_RR_dic, on="Vcf_id")
    plot_df = plot_df[plot_df['Cluster'] == 'C5']
    if domain in ['ABC', 'NB-ARC', 'AKR', 'P450', 'UGT']:
        if domain == 'NB-ARC':
            mask = (plot_df["Depth"] >= 0) & (plot_df["Depth"] <= 2)
        elif domain == 'P450':
            mask = (plot_df["Depth"] >= 0) & (plot_df["Depth"] <= 2.5)
        else:
            mask = (plot_df["Depth"] >= 0) & (plot_df["Depth"] <= 1.6)
        if domain == 'P450':
            sampled_in_range = plot_df[mask].sample(frac=0.02, random_state=42)  # 范围内采样
        else:
            sampled_in_range = plot_df[mask].sample(frac=0.05, random_state=42)  # 范围内采样
        out_of_range = plot_df[~mask]
        plot_df = pd.concat([sampled_in_range, out_of_range])

    sns.scatterplot(ax=ax_main, x='Als_GR50', y="Depth", linewidth=0, data=plot_df, color = Doamin_color[domain], alpha=0.5)
    sns.kdeplot(ax=ax_y, x=plot_df["Depth"], color = Doamin_color[domain], fill=True, vertical=True)
    if domain not in ['ABC', 'NB-ARC', 'AKR', 'P450']:
        ax_y.set_ylim([0, 2])
        ax_y.set_yticks(np.linspace(0, 2, 5))
        ax_main.set_ylim([0, 2])
        ax_main.set_yticks(np.linspace(0, 2, 5))
    ax_main.text(0.95, 0.95, domain,fontsize=15,
                     transform=ax_main.transAxes,
                     ha="right", va="top")
    ax_y.axhline(y=1,color='darkred',linestyle='--',linewidth=2)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['top'].set_visible(False)
    ax_main.set_xlabel("GR50 (g a.i. ha-1)", fontsize=12)
    ax_y.set_ylabel("")
    ax_y.spines['top'].set_visible(False)
    ax_y.spines['bottom'].set_visible(False)
    ax_y.spines['right'].set_visible(False)
    ax_y.xaxis.set_visible(False)
    ax_y.set_yticklabels([])


def draw_joyplot(depth_df, Vcf_RR_dic):
    depth_plot_df = pd.melt(
    depth_df,
    id_vars=['geneID'],
    var_name='Vcf_id',  # 新列名（存储转换前的列名）
    value_name='Depth'  # 新列名（存储转换前的值）
    )
    depth_plot_df['Type'] = depth_plot_df['Vcf_id'].map(Vcf_RR_dic)
    depth_plot_df.dropna(inplace=True)

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

def draw_boxplot(depth_df, Vcf_RR_dic):
    depth_plot_df = pd.melt(depth_df,id_vars=['geneID'],
                            var_name='Vcf_id',value_name='Depth')
    depth_plot_df['Type'] = depth_plot_df['Vcf_id'].map(Vcf_RR_dic)
    depth_plot_df.dropna(inplace=True)

    plt.figure(figsize=(3, 5))

    sns.boxplot(x="Type", y="Depth", data=depth_plot_df, palette=joy_color, showfliers=False)
    sns.stripplot(x="Type", y="Depth", data=depth_plot_df, palette=joy_color)

    plt.show()
    # # plt.savefig(os.path.join(out_dir,prefix+'_'+str(i*100)+'.pdf'), format='pdf')
    # plt.close()

def main(depth_file, herbicide_resistance_file, sample_list):

    fig = plt.figure(figsize=(10,9))
    outer = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    for idx, domain in enumerate(Doamin_color):
        cur_depth_file = os.path.join(depth_file, f"{domain}_uniq_ratio.txt")
        row, col = divmod(idx, 3)
        inner = outer[row, col].subgridspec(1,2, width_ratios=(6,1), hspace=0.02)
        ax_main = fig.add_subplot(inner[0,0])
        ax_y = fig.add_subplot(inner[0,1])

        depth_df = pd.read_table(cur_depth_file)
        depth_df.columns = ["geneID"] + depth_df.columns[1:].tolist()
        print(depth_df)

        HB_RR_Df = pd.read_table(herbicide_resistance_file, sep="\t").set_index("Vcf_id")
        Vcf_RR_dic = HB_RR_Df['Type'].to_dict()

        HB_RR_Df = HB_RR_Df.loc[sample_list, ]

        # draw_joyplot(depth_df, Vcf_RR_dic)
        # draw_boxplot(depth_df, Vcf_RR_dic)
        # draw_scatter(depth_df, HB_RR_Df, domain, ax)
    #     plot_custom_joint(ax_main, ax_y, depth_df, HB_RR_Df, domain)

    # plt.tight_layout()
    # # plt.show()
    # plt.savefig(rf'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\Domain_CN_domain_joint.pdf', format='pdf')

if __name__ == "__main__":
    # depth_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\P450_uniq_ratio.txt'
    herbicide_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
    filtered_sample_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\rmsus_sample.list'
    depth_file_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis'

    cur_sample_list = []
    with open(filtered_sample_file, 'r') as samplef:
        for i in samplef:
            line=i.strip()
            cur_sample_list.append(line)
    cur_sample_list.remove('15-18_1')
    cur_sample_list.remove('W-9_1')

    main(depth_file_dir, herbicide_resistance_file, cur_sample_list)