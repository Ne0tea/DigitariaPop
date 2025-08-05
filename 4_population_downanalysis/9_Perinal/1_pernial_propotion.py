'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-23 23:25:18
LastEditors: Ne0tea
LastEditTime: 2025-06-18 21:35:43
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'

colors = ["#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226"]
columns = [0, 1, 2, 3, 4,5]

def process_group(df, group_ids):
    group_df = df[df['ID'].isin(group_ids)].copy()
    sort_df = group_df.sort_values(columns)
    return sort_df

def draw_propotion(group1_df, group2_df, prefix):
    fig = plt.figure(figsize=(10, 4))
    gs = plt.GridSpec(2, 2, width_ratios=[4, 1])  # 2行2列，第1列宽度是第2列的4倍

    ax1 = fig.add_subplot(gs[0, 0])  # 第1行第1列
    ax3 = fig.add_subplot(gs[0, 1])  # 第1行第2列
    ax2 = fig.add_subplot(gs[1, 0])  # 第2行第1列
    ax4 = fig.add_subplot(gs[1, 1])  # 第2行第2列
    plt.subplots_adjust(hspace=0.4)
    # 绘制上方子图（组1）
    x1 = np.arange(len(group1_df))
    bottom1 = np.zeros(len(group1_df))
    for i, col in enumerate(columns):
        ax1.bar(x1, group1_df[col], bottom=bottom1, color=colors[i], label=col, width = 0.9)
        bottom1 += group1_df[col].values
    # ax1.set_title('Group 1 (Sorted by Total Value)', pad=20)
    ax1.set_ylabel('Value')
    ax1.set_xticks([])
    # ax1.set_xticklabels(group1_df['ID'], rotation=45, ha='right')

    x2 = np.arange(len(group2_df))
    bottom2 = np.zeros(len(group2_df))
    for i, col in enumerate(columns):
        ax2.bar(x2, group2_df[col], bottom=bottom2, color=colors[i], label=col, width = 0.9)
        bottom2 += group2_df[col].values
    # ax2.set_title('Group 2 (Sorted by Total Value)', pad=20)
    ax2.set_ylabel('Admixture proportion')
    ax2.set_xticks([])
    # ax2.set_xticklabels(group2_df['ID'], rotation=45, ha='right')

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, title='Categories',
               bbox_to_anchor=(1.05, 0.8), loc='upper left')

    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.set_yticks([])

    total_values = group1_df[columns].sum().values
    ax3.pie(total_values, colors=colors, labels=columns, startangle=90)

    total_values = group2_df[columns].sum().values
    ax4.pie(total_values, colors=colors, labels=columns, startangle=90)

    plt.tight_layout()
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_admixture'+r'\Pernial_propotion_'+prefix+'.pdf')

def main(proportion_file, non_499_sp_file, perinal_material_file):
    proportion_df = pd.read_csv(proportion_file, sep='\s+', header=None)
    non_499_sp_df = pd.read_csv(non_499_sp_file, sep='\t', header=None, names=['ID'])
    proportion_df['ID'] = non_499_sp_df['ID']

    perinal_ID_dic = {}
    with open(perinal_material_file, 'r') as f:
        for i in f:
            if i.startswith('ID'): continue
            line = i.strip().split('\t')
            perinal_ID = line[3]
            cur_year = str(line[4])
            material_ID = line[0]
            perinal_year = '20' + str(int(perinal_ID.split('_')[2]))
            if perinal_ID not in perinal_ID_dic:
                if cur_year == perinal_year:
                    perinal_ID_dic[perinal_ID] = {'Perinal':[material_ID]}
                else:
                    perinal_ID_dic[perinal_ID] = {'Current':[material_ID]}
            else:
                if cur_year == perinal_year:
                    if 'Perinal' not in perinal_ID_dic[perinal_ID]:
                        perinal_ID_dic[perinal_ID]['Perinal'] = [material_ID]
                    else:
                        perinal_ID_dic[perinal_ID]['Perinal'].append(material_ID)
                else:
                    if 'Current' not in perinal_ID_dic[perinal_ID]:
                        perinal_ID_dic[perinal_ID]['Current'] = [material_ID]
                    else:
                        perinal_ID_dic[perinal_ID]['Current'].append(material_ID)
    # print(pd.DataFrame(perinal_ID_dic).T.map(lambda x: len(x)))
    for i in perinal_ID_dic:
        print(i)
        perinal_ID_list =  perinal_ID_dic[i]['Perinal']
        current_ID_list = perinal_ID_dic[i]['Current']
        cur_perinal_df = process_group(proportion_df, perinal_ID_list)
        cur_current_df = process_group(proportion_df, current_ID_list)
        total_values = cur_perinal_df[columns].sum().apply(lambda x: x/len(perinal_ID_list))
        print(total_values)
        total_values = cur_current_df[columns].sum().apply(lambda x: x/len(current_ID_list))
        print(total_values)
    #     draw_propotion(cur_perinal_df, cur_current_df, i )

if __name__ == "__main__":
    proportion_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\P_structure\GATK3_nonsyn_499\Run2\GATK3_non_499_2.6.meanQ'
    non_499_sp_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\P_structure\GATK3_499.txt'
    perinal_material_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_group_material_loc.txt'
    main(proportion_file, non_499_sp_file, perinal_material_file)