'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-23 23:25:18
LastEditors: Ne0tea
LastEditTime: 2025-04-30 15:36:30
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'

# colors = ["#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226"]
colors = ["#005f73","#0a9396","#94d2bd","#001219","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226"]
columns = [0, 1, 2, 3, 4, 5]

def process_group(df, group_ids, perinal_ID_dic):
    group_df = df[df['ID'].isin(group_ids)].copy()
    sort_df = group_df.sort_values(columns)
    return sort_df

def draw_propotion(all_df, tick_loc, perinal_ID_dic):
    fig = plt.figure(figsize=(15, 2))
    gs = plt.GridSpec(1, len(tick_loc), width_ratios=[len(perinal_ID_dic[x]) for x in ['C4', 'Group_1', 'Group_2',\
        'Group_3', 'Group_4', 'Group_6','Group_7','Group_8','Group_12']], 
                      wspace=0.05)

    axC4 = fig.add_subplot(gs[0])
    groupC4_df = all_df.loc[perinal_ID_dic['C4']]
    ax1 = fig.add_subplot(gs[1])
    group1_df = all_df.loc[perinal_ID_dic['Group_1']]
    ax2 = fig.add_subplot(gs[2])
    group2_df = all_df.loc[perinal_ID_dic['Group_2']]
    ax3 = fig.add_subplot(gs[3])
    group3_df = all_df.loc[perinal_ID_dic['Group_3']]
    ax4 = fig.add_subplot(gs[4])
    group4_df = all_df.loc[perinal_ID_dic['Group_4']]
    # ax5 = fig.add_subplot(gs[5])
    # group5_df = all_df.loc[perinal_ID_dic['Group_5']]
    ax6 = fig.add_subplot(gs[5])
    group6_df = all_df.loc[perinal_ID_dic['Group_6']]
    ax7 = fig.add_subplot(gs[6])
    group7_df = all_df.loc[perinal_ID_dic['Group_7']]
    ax8 = fig.add_subplot(gs[7])
    group8_df = all_df.loc[perinal_ID_dic['Group_8']]
    ax12 = fig.add_subplot(gs[8])
    group12_df = all_df.loc[perinal_ID_dic['Group_12']]

    x1 = np.arange(len(group1_df))
    bottom1 = np.zeros(len(group1_df))
    for i, col in enumerate(columns):
        ax1.bar(x1, group1_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group1_df[col].values
    ax1.set_title('Group_1')

    x1 = np.arange(len(group2_df))
    bottom1 = np.zeros(len(group2_df))
    for i, col in enumerate(columns):
        ax2.bar(x1, group2_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group2_df[col].values
    ax2.set_title('Group_2')

    x1 = np.arange(len(group3_df))
    bottom1 = np.zeros(len(group3_df))
    for i, col in enumerate(columns):
        ax3.bar(x1, group3_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group3_df[col].values
    ax3.set_title('Group_3')

    x1 = np.arange(len(group4_df))
    bottom1 = np.zeros(len(group4_df))
    for i, col in enumerate(columns):
        ax4.bar(x1, group4_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group4_df[col].values
    ax4.set_title('Group_4')

    # x1 = np.arange(len(group5_df))
    # bottom1 = np.zeros(len(group5_df))
    # for i, col in enumerate(columns):
    #     ax5.bar(x1, group5_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
    #     bottom1 += group5_df[col].values
    # ax5.set_title('Group_5')

    x1 = np.arange(len(group6_df))
    bottom1 = np.zeros(len(group6_df))
    for i, col in enumerate(columns):
        ax6.bar(x1, group6_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group6_df[col].values
    ax6.set_title('Group_6')

    x1 = np.arange(len(group7_df))
    bottom1 = np.zeros(len(group7_df))
    for i, col in enumerate(columns):
        ax7.bar(x1, group7_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group7_df[col].values
    ax7.set_title('Group_7')

    x1 = np.arange(len(group8_df))
    bottom1 = np.zeros(len(group8_df))
    for i, col in enumerate(columns):
        ax8.bar(x1, group8_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group8_df[col].values
    ax8.set_title('Group_8')
    print(group1_df)
    # ax8.set_xticklabels(group8_df.index, rotation=90, ha='right')

    x1 = np.arange(len(group12_df))
    bottom1 = np.zeros(len(group12_df))
    for i, col in enumerate(columns):
        ax12.bar(x1, group12_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += group12_df[col].values
    ax12.set_title('Group_12')

    x1 = np.arange(len(groupC4_df))
    bottom1 = np.zeros(len(groupC4_df))
    for i, col in enumerate(columns):
        axC4.bar(x1, groupC4_df[col], bottom=bottom1, color=colors[i], label=col, width = 1)
        bottom1 += groupC4_df[col].values
    axC4.set_title('C4')

    for ax in [ax1, ax2,ax3, ax4,ax6,ax7, ax8,ax12, axC4]:
        ax.xaxis.set_visible(False)
        if ax != axC4:
            ax.yaxis.set_visible(False)
            ax.spines['left'].set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.tight_layout()
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Sympatric_propotion_all_p6.pdf')

def main(proportion_file, non_499_sp_file, perinal_material_file):
    proportion_df = pd.read_csv(proportion_file, sep='\s+', header=None)
    non_499_sp_df = pd.read_csv(non_499_sp_file, sep='\t', header=None, names=['ID'])
    proportion_df['ID'] = non_499_sp_df['ID']
    proportion_df.set_index('ID', inplace=True)
    print(proportion_df)

    perinal_ID_dic = {}
    with open(perinal_material_file, 'r') as f:
        for i in f:
            if i.startswith('Year'): continue
            line = i.strip().split('\t')
            perinal_ID = line[7]
            material_ID = line[1]
            if line[6] == 'C4': 
                perinal_ID = 'C4'
            if perinal_ID not in perinal_ID_dic:
                perinal_ID_dic[perinal_ID] = [material_ID]
            else:
                perinal_ID_dic[perinal_ID].append(material_ID)

    tick_loc = []
    all_material = []
    left = 0
    for i in perinal_ID_dic:
        cur_material = perinal_ID_dic[i]
        cur_adm_df = proportion_df.loc[cur_material]
        # print(len(cur_adm_df) // 2)
        tick_loc.append(len(cur_adm_df) // 2 + left)
        left += len(cur_adm_df)
        all_material.extend(cur_material)

    all_material_df = proportion_df.loc[all_material]
    draw_propotion(all_material_df, tick_loc, perinal_ID_dic)

if __name__ == "__main__":
    # proportion_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\P_structure\GATK3_nonsyn_499\Run2\GATK3_non_499_2.6.meanQ'
    proportion_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\P_structure\GATK3\faststructure\GATK3.6.meanQ'
    non_499_sp_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\P_structure\Dsan_all_sample.list'
    sympatric_material_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Sympatric_material2.txt'
    main(proportion_file, non_499_sp_file, sympatric_material_file)