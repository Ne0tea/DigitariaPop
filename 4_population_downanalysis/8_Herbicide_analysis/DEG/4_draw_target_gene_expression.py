'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-17 20:45:29
LastEditors: Ne0tea
LastEditTime: 2025-08-01 12:18:28
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from matplotlib.ticker import FixedLocator
import re

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
treat_color={"S": "#f2c22f", "R": "#2d674d"}
custom_orders = {
    'C': ['0h', '3h', '12h'],
    'X': ['0h', '24h', '48h'],
    'Y': ['0h', '6h', '24h'],
}
compare_pairs = {'C':[((0, 'R'), (1, 'R')), ((0, 'R'), (2, 'R')),
                      ((0, 'S'), (1, 'S')), ((0, 'S'), (2, 'S'))],
                 'X':[((0, 'R'), (1, 'R')), ((0, 'R'), (2, 'R')),
                      ((0, 'S'), (1, 'S')), ((0, 'S'), (2, 'S'))],
                 'Y':[((0, 'R'), (1, 'R')), ((0, 'R'), (2, 'R')),
                      ((0, 'S'), (1, 'S')), ((0, 'S'), (2, 'S')),
                      ((0, 'R'), (0, 'S')), ((1, 'R'), (1, 'S')),((2, 'R'), (2, 'S'))]
}
def draw_bar_plot(detailed_data, target_gene_name):
    time_ord = detailed_data['Time'].sort_values().unique()
    time_ord = sorted(time_ord, key=lambda x: int(re.search(r'\d+', x).group()))
    treat_ord = ['C','X','Y']
    expr_plot = sns.catplot(x="day_ordered",y="FPKM",hue="Var",row="Ht",
                            row_order=treat_ord,
                            data=detailed_data,
                            kind="bar",
                            aspect=1,
                            height=3,
                            linewidth=1.5,
                            palette=treat_color,
                            errorbar=("pi", 95),
                            capsize=0.5,
                            fill=False,
                            err_kws={'linewidth': 1},
                            # width=0.4,
                            # order=time_ord,
                            sharex=False,
                            sharey=False
                        )
    plt.subplots_adjust(
        hspace=0.05   # 垂直间距
    )

    for ax, treat in zip(expr_plot.axes.flat, treat_ord):
        locs = ax.get_xticks()
        ax.xaxis.set_major_locator(FixedLocator(locs))
        sns.stripplot(x="day_ordered",y="FPKM",hue="Var",data=detailed_data[detailed_data['Ht'] == treat],
                      dodge=True,palette=treat_color,alpha=0.5,ax=ax)
        ax.set_xticklabels(custom_orders[treat])
        annotator = Annotator(ax, compare_pairs[treat], data=detailed_data[detailed_data['Ht'] == treat], x="day_ordered",y="FPKM",
                                hue="Var")
        annotator.configure(test="t-test_ind",text_format="star",loc="inside",
                            # hide_non_significant=True,
                            # comparisons_correction="benjamini-hochberg"
                            )
        annotator.apply_and_annotate()
        if ax.legend_:
            ax.legend_.remove()
    expr_plot.set(xlabel="")
    expr_plot.set_titles(row_template="{row_name}",size=10,color="black")
    plt.tight_layout()
    # plt.savefig(f'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\{target_gene_name}_DE.pdf')
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
            # if srr_id == '15-2-0-3': continue
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

def main(treat_file, target_gene_name, R_count_file, S_count_file):
    R_count_df = read_stringtie_out(R_count_file)
    S_count_df = read_stringtie_out(S_count_file)
    # S_count_df = S_count_df.drop(columns=['15-2-0-3'])

    treat_name_SRR, CK_list = make_treat_dic(treat_file)
    # CK_list.remove('15-2-0-3')

    R_count_df = R_count_df.loc[[target_gene_name]]
    S_count_df = S_count_df.loc[[target_gene_name]]

    all_treat_name = []
    for i in treat_name_SRR:
        all_treat_name.extend(treat_name_SRR[i])

    R_long_df = get_treatment_expr(R_count_df, all_treat_name, CK_list, '21-17')
    S_long_df = get_treatment_expr(S_count_df, all_treat_name, CK_list, '15-2')
    All_long_df = pd.concat([R_long_df, S_long_df], axis=0).reset_index(drop=True)
    All_long_df['Var'] = All_long_df['Var'].replace('21', 'R')
    All_long_df['Var'] = All_long_df['Var'].replace('15', 'S')

    CK_treat_row = (All_long_df['Time'] == '0h') & (All_long_df['Ht'] == 'C')
    CK_treat_df = All_long_df[CK_treat_row].copy()
    CK_treat_X_df = CK_treat_df.assign(Ht='X')
    CK_treat_Y_df = CK_treat_df.assign(Ht='Y')
    All_long_df = pd.concat([All_long_df, CK_treat_X_df, CK_treat_Y_df], ignore_index=True)

    All_long_df['day_ordered'] = All_long_df.apply(
        lambda row: custom_orders[row['Ht']].index(row['Time'])
                    if row['Time'] in custom_orders[row['Ht']]
                    else len(custom_orders[row['Time']]),
        axis=1
    )
    draw_bar_plot(All_long_df, target_gene_name)

if __name__ == "__main__":
    treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat.txt'
    target_gene_name = 'Chr13.1090'
    R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_fpkm_matrix.csv'
    S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_fpkm_matrix.csv'

    main(treat_file, target_gene_name, R_count_file, S_count_file)
