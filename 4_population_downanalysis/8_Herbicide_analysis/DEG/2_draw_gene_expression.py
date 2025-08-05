'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-17 20:45:29
LastEditors: Ne0tea
LastEditTime: 2025-08-01 15:27:30
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import re

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
treat_color={"R": "#2d674d", "S": "#f2c22f"}

compare_pairs = {'C':[(('0h', 'R'), ('6h', 'R')), (('0h', 'R'), ('24h', 'R')),
                      (('0h', 'S'), ('6h', 'S')), (('0h', 'S'), ('24h', 'S'))],
                 'X':[(('0h', 'R'), ('6h', 'R')), (('0h', 'R'), ('24h', 'R')),
                      (('0h', 'S'), ('6h', 'S')), (('0h', 'S'), ('24h', 'S'))],
                 'Y':[(('0h', 'R'), ('6h', 'R')), (('0h', 'R'), ('24h', 'R')),
                      (('0h', 'S'), ('6h', 'S')), (('0h', 'S'), ('24h', 'S')),
                      (('0h', 'R'), ('0h', 'S')), (('6h', 'R'), ('6h', 'S')),(('24h', 'R'), ('24h', 'S'))]
}

def draw_bar_plot(detailed_data, gene_name_ord, treatname):
    time_ord = detailed_data['Time'].sort_values().unique()
    time_ord = sorted(time_ord, key=lambda x: int(re.search(r'\d+', x).group()))

    expr_plot = sns.catplot(x="Time",y="FPKM",hue="Var",row="gene_id",
                            row_order=gene_name_ord,
                            data=detailed_data,
                            kind="bar",
                            aspect=3,
                            height=1,
                            palette=treat_color,
                            errorbar=("pi", 95),
                            capsize=0.5,
                            fill=False,
                            err_kws={'linewidth': 1},
                            # width=0.6,
                            sharey=False,
                            order=time_ord
                        )
    for ax, gene_id in zip(expr_plot.axes.flat, gene_name_ord):
        sns.swarmplot(x="Time",y="FPKM",hue="Var",data=detailed_data[detailed_data['gene_id'] == gene_id],dodge=True,palette=treat_color,alpha=0.5,
            ax=ax
        )
        annotator = Annotator(ax, compare_pairs[treatname], data=detailed_data[detailed_data['Ht'] == treatname], x="Time",y="FPKM",
                                hue="Var")
        annotator.configure(test="t-test_ind",text_format="star",loc="inside",
                            hide_non_significant=True,
                            # comparisons_correction="benjamini-hochberg"
                            )
        annotator.apply_and_annotate()
        if ax.legend_:
            ax.legend_.remove()

    expr_plot.set_titles(row_template="{row_name}"+'_'+treatname,size=10,color="black")
    plt.tight_layout()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\DsALS' + '_' + treatname+'.pdf')
    plt.show()

def read_stringtie_out(file):
    result_df = pd.read_csv(file, header=0)
    result_df = result_df.set_index('gene_id')
    return result_df

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

def get_treatment_expr(raw_expr_df, all_treat_name, CK_list, prefix):
    temp_Sen_id = [x for x in all_treat_name if x.startswith(prefix)]

    cur_count_df = raw_expr_df[temp_Sen_id]
    cur_count_df = cur_count_df.reset_index()
    cur_count_df_long = cur_count_df.melt(id_vars="gene_id",var_name="sample",value_name="FPKM")
    cur_count_df_long[["Var", "Ht", "Time", "Rep"]] = cur_count_df_long["sample"].str.split("-", expand=True)[[0, 2,3,4]]

    cur_CK_df = raw_expr_df[[x for x in CK_list if x.startswith(prefix)]]
    cur_CK_df = cur_CK_df.reset_index()
    cur_CK_df_long = cur_CK_df.melt(id_vars="gene_id",var_name="sample",value_name="FPKM")
    cur_CK_df_long[["Var", "Time", "Rep"]] = cur_CK_df_long["sample"].str.split("-", expand=True)[[0, 2,3]]
    cur_CK_df_long['Ht'] = cur_count_df_long["Ht"].to_list()[0]
    cur_CK_df_long['Time'] = '0h'
    cur_count_df_long = pd.concat([cur_count_df_long, cur_CK_df_long], axis=0)

    return cur_count_df_long

def main(treat_file, target_gene_file, R_count_file, S_count_file):
    R_count_df = read_stringtie_out(R_count_file)
    S_count_df = read_stringtie_out(S_count_file)
    # S_count_df = S_count_df.drop(columns=['15-2-0-3'])

    treat_name_SRR, CK_list = make_treat_dic(treat_file)
    # CK_list.remove('15-2-0-3')

    target_gene_list = []
    with open(target_gene_file, 'r') as tf:
        for i in tf:
            line=i.strip().split()
            target_gene_list.append(line[0])
    R_count_df = R_count_df.loc[target_gene_list]
    S_count_df = S_count_df.loc[target_gene_list]

    for treatname in treat_name_SRR:
        if treatname != 'Y': continue
        cur_R_long_df = get_treatment_expr(R_count_df, treat_name_SRR[treatname], CK_list, '21-17')
        cur_S_long_df = get_treatment_expr(S_count_df, treat_name_SRR[treatname], CK_list, '15-2')
        cur_all_df = pd.concat([cur_R_long_df, cur_S_long_df], axis=0).reset_index(drop=True)

        cur_all_df['Var'] = cur_all_df['Var'].replace('21', 'R')
        cur_all_df['Var'] = cur_all_df['Var'].replace('15', 'S')
        max_fpkm = cur_all_df.groupby('gene_id')['FPKM'].max()
        # cur_all_df = cur_all_df.groupby('gene_id').filter(lambda x: x['FPKM'].max() > 10)

        low_expr_gene = max_fpkm[max_fpkm < 5].index.tolist()
        high_expr_gene = max_fpkm[max_fpkm >= 5].index.tolist()
        strings = '\t'.join(low_expr_gene)
        print(f'Current low expr gene id {strings}')
        # print(cur_all_df)
        draw_bar_plot(cur_all_df, high_expr_gene, treatname)

if __name__ == "__main__":
    treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat.txt'
    target_gene_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_target_gene.list'
    R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_fpkm_matrix.csv'
    S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_fpkm_matrix.csv'

    main(treat_file, target_gene_file, R_count_file, S_count_file)
