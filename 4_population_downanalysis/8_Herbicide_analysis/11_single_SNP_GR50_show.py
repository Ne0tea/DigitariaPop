'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-07 22:53:16
LastEditors: Ne0tea
LastEditTime: 2025-07-13 17:05:01
'''
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import re

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
outlier=['MMT-13-20_1']
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
hap_color={'0':'#f5b940','1':'#b6d6d6','2':'#005d5c'}
def custom_replace(column_values, replacement_dict):
    return column_values.replace(replacement_dict)
def hap_key(x):
    return int(re.search(r'Hap(\d+)', x).group(1))
def draw_plot_boxplot(plot_df, prefix):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    plot_df = plot_df[plot_df['Target_loc'] != -1]
    plot_df['Target_loc'] = plot_df['Target_loc'].astype(str)
    plt.figure(figsize=(4, 2))
    box=sns.boxplot(
        x="Target_loc",
        y="Herbicide",
        data=plot_df,
        width=0.5,
        # orient = 'h',
        hue='Target_loc',
        fill=False,
        showfliers=False,
        palette=hap_color
    )
    anno = Annotator(ax = box,data=plot_df, x='Target_loc', y='Herbicide',
                     pairs=[('2','0'), ('1','0')], 
                     permutations=1000,
                     orient='v', verbose=False)
    anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
    anno.apply_and_annotate()
    sns.stripplot(x="Target_loc", y='Herbicide', data=plot_df,
                  palette=hap_color,
                    legend=False,
                    dodge=False, jitter=True, alpha=0.5, edgecolor='black')
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    # plt.show()
    plt.savefig(prefix+'_GR50_boxplot.pdf', format='pdf')

def draw_Ecotype_propotion(plot_df, prefix):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    target_column = 'Resistance'
    plot_df = plot_df[plot_df['Target_loc'] != -1]
    plot_df['Target_loc'] = plot_df['Target_loc'].astype(str)
    plot_df = plot_df.groupby(['Target_loc', target_column]).size().reset_index(name='Count')
    plot_df = plot_df.pivot(index='Target_loc', columns=target_column, values='Count')
    plot_df.fillna(0, inplace=True)

    row_sums = plot_df.sum(axis=0)
    plot_df = (plot_df.div(row_sums, axis=1) * 100).T

    Clades = plot_df.index.tolist()
    haps = plot_df.columns.tolist()
    bottom = [0] * len(Clades)

    fig, ax = plt.subplots(figsize=(4,2))
    for i, hap in enumerate(haps):
        counts = plot_df[hap].values
        ax.barh(Clades, counts, left=bottom, label=hap, color=hap_color[hap])
        bottom = [b + c for b, c in zip(bottom, counts)]

    for i, (name, value) in enumerate(row_sums.items()):
        ax.text(1.02, i, f"N={int(value)}", 
                va='center', ha='left', 
                transform=ax.get_yaxis_transform())

    plt.title("")
    plt.xlim(-5, 105)
    plt.ylabel("")
    plt.xlabel("Percentage (%)")
    # plt.legend(ncol=len(Clades)/2)
    # plt.show()
    plt.savefig(prefix+'_Eco_Prop.pdf', format='pdf')

def main(hap_file, indv_file, target_snp_loc, pos_file, GR50_sample_file, sample_clade_file, sus_sample_file):
    sus_sample_list = []
    with open(sus_sample_file, 'r') as sus_sample_f:
        for i in sus_sample_f:
            line=i.strip().split()
            sus_sample_list.append(line[0])

    indv_list = []
    with open(indv_file, 'r') as indv_f:
        for i in indv_f:
            line=i.strip().split()
            indv_list.append(line[0])

    pos_list = []
    with open(pos_file, 'r') as pos_f:
        for i in pos_f:
            line=i.strip().split('\t')
            pos_list.append(line[1])

    hap_012_df = pd.read_table(hap_file, header=None, index_col=0)
    hap_012_df.index = indv_list
    hap_012_df.columns = pos_list
    hap_target_loc = hap_012_df[target_snp_loc]

    HR_mat_df = pd.read_table(GR50_sample_file)
    HR_mat_df.index = indv_list
    Clade_mat_df = pd.read_table(sample_clade_file, header=0)
    Clade_mat_df.index = indv_list
    Clade_dic = Clade_mat_df['Resistance'].to_dict()

    hap_HR_df = pd.concat([hap_target_loc, HR_mat_df, Clade_mat_df], axis=1)
    # hap_HR_df.dropna(axis=0, inplace=True)
    hap_HR_df = hap_HR_df[hap_HR_df['Clade'].isin(['C5'])]
    # hap_HR_df = hap_HR_df[hap_HR_df['Clade'].isin(['C5']) & ~hap_HR_df['sample'].isin(sus_sample_list)]
    hap_HR_df.dropna(inplace=True)
    hap_HR_df = hap_HR_df.rename(columns={target_snp_loc:'Target_loc'})

    print(hap_HR_df.groupby(['Target_loc']).count())
    # draw_plot_boxplot(hap_HR_df, hap_file)
    # draw_Ecotype_propotion(hap_HR_df, hap_file)
    # for col in target_snp_hap_df.columns:
    #     target_snp_hap_df[col] = custom_replace(target_snp_hap_df[col], col_replace_dic.get(str(col), {}))
    # target_snp_hap_df.columns = [f"{col}_{snp_eff_dic.get(str(col), {})[1]}" for col in target_snp_hap_df.columns]
    # target_snp_hap_df.to_csv(hap_file+'_hap.txt', sep='\t')

if __name__ == '__main__':
    # snp_var_file = sys.argv[1]
    # hap_file = sys.argv[2]
    # indv_file = sys.argv[3]
    # pos_file = sys.argv[4]
    # GR50_sample_file = sys.argv[5]
    # sample_clade_file = sys.argv[6]

    # snp_var_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Chr20_1762_region.1534.gene.snp"
    hap_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Chr20_1762_case\Chr20_temp.012"
    indv_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Chr20_1762_case\Chr20_temp.012.indv"
    pos_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Chr20_1762_case\Chr20_temp.012.pos"
    GR50_sample_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_Herbicide_GR50.txt"
    # sample_clade_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_clade.txt'
    sample_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Sample_Resistance_level.list'
    sus_sample_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Suspective_sample_nearGR5090.list'
    target_snp_loc = '34032644'
    main(hap_file, indv_file, target_snp_loc, pos_file, GR50_sample_file, sample_resistance_file, sus_sample_file)
