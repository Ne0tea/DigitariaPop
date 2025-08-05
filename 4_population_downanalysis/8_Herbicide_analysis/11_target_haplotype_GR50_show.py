'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-07 22:53:16
LastEditors: Ne0tea
LastEditTime: 2025-05-27 14:50:55
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
hap_color={'Hap1':'#f5b940','Hap2':'#005d5c'}
def custom_replace(column_values, replacement_dict):
    return column_values.replace(replacement_dict)
def hap_key(x):
    return int(re.search(r'Hap(\d+)', x).group(1))
def draw_plot_boxplot(plot_df, prefix):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    plot_df = plot_df[plot_df['Clade'] == 'C5']
    plot_df = plot_df[~plot_df.index.isin(outlier)]
    plot_df = plot_df[plot_df['Hap'].isin(['Hap1','Hap2'])]
    plt.figure(figsize=(2, 3))
    box=sns.boxplot(
        x="Hap",
        y="Herbicide",
        data=plot_df,
        width=0.5,
        # orient = 'h',
        hue='Hap',
        fill=False,
        showfliers=False,
        palette=hap_color
    )
    anno = Annotator(ax = box,data=plot_df, x='Hap', y='Herbicide',
                     pairs=[('Hap1','Hap2')], 
                     permutations=1000,
                     orient='v', verbose=False)
    anno.configure(test='t-test_ind', text_format='simple', line_height=0.03, line_width=1, hide_non_significant=True)
    anno.apply_and_annotate()
    sns.stripplot(x="Hap", y='Herbicide', data=plot_df,
                  palette=hap_color,
                    legend=False,
                    dodge=False, jitter=True, alpha=0.5, edgecolor='black')
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
    # plt.savefig(prefix+'_GR50_boxplot.pdf', format='pdf')

def draw_Ecotype_propotion(plot_df, prefix):
    print(plot_df)
    target_column = 'Resistance'
    plot_df.loc['MMT-13-20_1', target_column] = 'R_2'
    plot_df.loc['MMT-13-10_1', target_column] = 'R_2'
    plot_df = plot_df[plot_df['Hap'].isin(['Hap1','Hap2'])]

    plot_df = plot_df.groupby(['Hap', target_column]).size().reset_index(name='Count')
    plot_df = plot_df.pivot(index='Hap', columns=target_column, values='Count')
    plot_df.fillna(0, inplace=True)

    sorted_index = sorted(plot_df.index, key=hap_key, reverse=True)
    plot_df = plot_df.reindex(index=sorted_index).T
    row_sums = plot_df.sum(axis=1)
    plot_df = (plot_df.div(row_sums, axis=0) * 100).T
    print(plot_df)

    Clades = plot_df.columns.tolist()
    haps = plot_df.index.tolist()
    bottom = [0] * len(Clades)

    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    fig, ax = plt.subplots(figsize=(5,2))
    for i, hap in enumerate(haps):
        counts = plot_df.loc[hap].values
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
    plt.show()
    # plt.savefig(prefix+'_Eco_Prop.pdf', format='pdf')

def main(hap_file, indv_file, pos_file, GR50_sample_file, sample_clade_file, sus_sample_file):
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
    hap_012_df = hap_012_df[[col for col in hap_012_df.columns if (hap_012_df[col] == 2).sum() != 1]]

    HR_mat_df = pd.read_table(GR50_sample_file)
    HR_mat_df.index = indv_list
    Clade_mat_df = pd.read_table(sample_clade_file, header=0)
    Clade_mat_df.index = indv_list
    Clade_dic = Clade_mat_df['Resistance'].to_dict()

    hap_HR_df = pd.concat([hap_012_df, HR_mat_df, Clade_mat_df], axis=1)
    # hap_HR_df.dropna(axis=0, inplace=True)
    hap_HR_df = hap_HR_df[hap_HR_df['Clade'].isin(['C5']) & ~hap_HR_df['sample'].isin(sus_sample_list)]
    hap_HR_df.dropna(inplace=True)

    ## filtered snp and material
    target_snp_df = hap_HR_df.iloc[:,:-5]
    # target_snp_df.replace(to_replace=[-1,1], value=[0,0], inplace=True)
    target_snp_hap_df = target_snp_df.drop_duplicates(keep='first')
    target_snp_hap_df = target_snp_hap_df.reset_index(drop=True)
    target_snp_hap_df.index = ['Hap' + str(i+1) for i in range(len(target_snp_hap_df))]
    print(target_snp_hap_df)

    hap_GR50_df = pd.DataFrame()
    hap_Eco_df = pd.DataFrame()
    for row,row_df in target_snp_hap_df.iterrows():
        matched_hap_rows = target_snp_df[target_snp_df.eq(row_df).all(axis=1)]
        matched_hap_rows_index = matched_hap_rows.index.to_list()
        matched_hap_GR50 = hap_HR_df[hap_HR_df.index.isin(matched_hap_rows_index)]['Herbicide'].to_frame()
        matched_hap_GR50['Hap'] = row
        matched_hap_GR50['Clade'] = matched_hap_GR50.index.map(Clade_dic)
        hap_GR50_df = pd.concat([hap_GR50_df, matched_hap_GR50], axis=0)

        matched_hap_Eco = hap_HR_df[hap_HR_df.index.isin(matched_hap_rows_index)]['Resistance'].to_frame()
        matched_hap_Eco['Hap'] = row
        matched_hap_Eco['Resistance'] = matched_hap_Eco.index.map(Clade_dic)
        hap_Eco_df = pd.concat([hap_Eco_df, matched_hap_Eco], axis=0)
    print(hap_Eco_df[(hap_Eco_df['Resistance'] == 'R_0') & (hap_Eco_df['Hap'] == 'Hap2')])
    print(hap_Eco_df.groupby('Resistance')['Hap'].value_counts())
    # print(hap_Eco_df)

    # draw_plot_boxplot(hap_GR50_df, hap_file)
    draw_Ecotype_propotion(hap_Eco_df, hap_file)
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
    hap_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr20.1762.mis.012"
    indv_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr20.1762.mis.012.indv"
    pos_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr20.1762.mis.012.pos"
    GR50_sample_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_Herbicide_GR50.txt"
    # sample_clade_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_clade.txt'
    sample_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Sample_Resistance_level.list'
    sus_sample_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Suspective_sample_nearGR5090.list'
    main(hap_file, indv_file, pos_file, GR50_sample_file, sample_resistance_file, sus_sample_file)
