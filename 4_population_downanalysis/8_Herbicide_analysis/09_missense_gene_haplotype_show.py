'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-07 22:53:16
LastEditors: Ne0tea
LastEditTime: 2025-04-10 16:49:46
'''
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
def custom_replace(column_values, replacement_dict):
    return column_values.replace(replacement_dict)
def hap_key(x):
    return int(re.search(r'Hap(\d+)', x).group(1))
def draw_plot_boxplot(plot_df, prefix):
    sns.set_theme(style="darkgrid")
    plt.figure(figsize=(6, 3))
    sns.boxplot(
        x="Herbicide",
        y="Hap",
        data=plot_df,
        width=0.7,
        orient = 'h',
        hue='Hap',
        fill=False,
        showfliers=False,
        palette='Blues'
        # palette="Set2" # 可选颜色
    )
    sns.stripplot(x='Herbicide', y="Hap", data=plot_df,
                  color='#003049',
                    legend=False,
                    dodge=True, jitter=True, alpha=0.5, edgecolor='black')
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    # plt.show()
    plt.savefig(prefix+'_GR50_boxplot.pdf', format='pdf')

def draw_Ecotype_propotion(plot_df, prefix):
    plot_df = plot_df.groupby(['Hap', 'Ecotype']).size().reset_index(name='Count')
    plot_df = plot_df.pivot(index='Hap', columns='Ecotype', values='Count')
    plot_df.fillna(0, inplace=True)

    sorted_index = sorted(plot_df.index, key=hap_key, reverse=True)
    plot_df = plot_df.reindex(index=sorted_index)
    row_sums = plot_df.sum(axis=1)
    plot_df = plot_df.div(row_sums, axis=0) * 100

    ecotypes = plot_df.columns.tolist()
    haps = plot_df.index.tolist()
    bottom = [0] * len(haps)

    sns.set_theme(style="darkgrid")
    fig, ax = plt.subplots(figsize=(10,6))
    for i, ecotype in enumerate(ecotypes):
        counts = plot_df[ecotype].values
        ax.barh(haps, counts, left=bottom, label=ecotype, color=Eco_type_color[ecotype])
        bottom = [b + c for b, c in zip(bottom, counts)]

    for i, (name, value) in enumerate(row_sums.items()):
        ax.text(1.02, i, f"N={int(value)}", 
                va='center', ha='left', 
                transform=ax.get_yaxis_transform())

    plt.title("")
    plt.xlim(-5, 105)
    plt.ylabel("")
    plt.xlabel("Percentage (%)")
    plt.legend(ncol=len(haps)/2)
    # plt.show()
    plt.savefig(prefix+'_Eco_Prop.pdf', format='pdf')

def main(snp_var_file, hap_file, indv_file, pos_file, GR50_sample_file, sample_clade_file):
    snp_eff_dic = {}
    with open(snp_var_file, 'r') as snp_var_f:
        for i in snp_var_f:
            line=i.strip().split('\t')
            fields = line[7].split('|')
            ref_snp, alt_snp = line[3], line[4]
            p_change = fields[10]
            if not p_change: continue
            his33arg = p_change.split('.')[1]
            alt_type = fields[1]
            snp_eff_dic[line[1]] = [(ref_snp, alt_snp), his33arg, alt_type]
    col_replace_dic = {x:{0:snp_eff_dic[x][0][0], 2:snp_eff_dic[x][0][1]} for x in snp_eff_dic}
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

    HR_mat_df = pd.read_table(GR50_sample_file)
    HR_mat_df.index = indv_list
    Clade_mat_df = pd.read_table(sample_clade_file)
    Clade_mat_df.index = indv_list

    hap_HR_df = pd.concat([hap_012_df, HR_mat_df, Clade_mat_df], axis=1)
    # hap_HR_df.dropna(axis=0, inplace=True)
    hap_HR_df = hap_HR_df[hap_HR_df['Clade']=='C5']

    target_snp = []
    for pos in pos_list:
        if pos not in snp_eff_dic: continue
        target_snp.append(pos)

    ## filtered snp and material
    target_snp_df = hap_HR_df[target_snp]
    target_snp_df.replace(to_replace=[-1,1], value=[0,0], inplace=True)
    target_snp_hap_df = target_snp_df.drop_duplicates(keep='first')
    target_snp_hap_df = target_snp_hap_df.reset_index(drop=True)
    target_snp_hap_df.index = ['Hap' + str(i+1) for i in range(len(target_snp_hap_df))]

    hap_GR50_df = pd.DataFrame()
    hap_Eco_df = pd.DataFrame()
    for row,row_df in target_snp_hap_df.iterrows():
        matched_hap_rows = target_snp_df[target_snp_df.eq(row_df).all(axis=1)]
        matched_hap_rows_index = matched_hap_rows.index.to_list()
        matched_hap_GR50 = hap_HR_df[hap_HR_df.index.isin(matched_hap_rows_index)]['Herbicide'].to_frame()
        matched_hap_GR50['Hap'] = row
        hap_GR50_df = pd.concat([hap_GR50_df, matched_hap_GR50], axis=0)

        matched_hap_Eco = hap_HR_df[hap_HR_df.index.isin(matched_hap_rows_index)]['Ecotype'].to_frame()
        matched_hap_Eco['Hap'] = row
        hap_Eco_df = pd.concat([hap_Eco_df, matched_hap_Eco], axis=0)

    draw_plot_boxplot(hap_GR50_df, hap_file)
    draw_Ecotype_propotion(hap_Eco_df, hap_file)
    for col in target_snp_hap_df.columns:
        target_snp_hap_df[col] = custom_replace(target_snp_hap_df[col], col_replace_dic.get(str(col), {}))
    target_snp_hap_df.columns = [f"{col}_{snp_eff_dic.get(str(col), {})[1]}" for col in target_snp_hap_df.columns]
    target_snp_hap_df.to_csv(hap_file+'_hap.txt', sep='\t')
if __name__ == '__main__':
    # snp_var_file = sys.argv[1]
    # hap_file = sys.argv[2]
    # indv_file = sys.argv[3]
    # pos_file = sys.argv[4]
    # GR50_sample_file = sys.argv[5]
    # sample_clade_file = sys.argv[6]

    snp_var_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Chr12.1534.gene.snp"
    hap_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr12.1534.012"
    indv_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr12.1534.012.indv"
    pos_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_candidate_gene\Chr12.1534.012.pos"
    GR50_sample_file = r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_Herbicide_GR50.txt"
    sample_clade_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\Vcf_sample_clade.txt'
    main(snp_var_file, hap_file, indv_file, pos_file, GR50_sample_file, sample_clade_file)
