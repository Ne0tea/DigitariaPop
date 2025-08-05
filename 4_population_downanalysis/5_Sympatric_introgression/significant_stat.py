'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-24 17:35:36
LastEditors: Ne0tea
LastEditTime: 2025-06-04 19:55:21
'''
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statannotations.Annotator import Annotator

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

def read_fd_df(fd_file):
    fd_file_name = os.path.basename(fd_file)
    run_name = re.search(r'^(RUN\d+)_', fd_file_name).group(1)
    group_name = re.search(r'_(Group_\d+)_', fd_file_name).group(1)
    print(group_name,run_name)

    if os.path.isfile(fd_file):
        cur_df=pd.read_csv(fd_file,header=0)
        cur_df = cur_df.drop_duplicates()
        cur_df['group']=group_name
        cur_df['run']=run_name
    return cur_df

Chr_subg_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\0_Dsan_chr_trans'
Chr_subg_df = pd.read_table(Chr_subg_file,names=['scaffold','Subg_Chr','Subg'])

sign_window_df=pd.DataFrame()

# target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Significant'
# for file_name in os.listdir(target_directory):
#     if file_name.endswith('10kwindows.txt'):
#         file_path = os.path.join(target_directory, file_name)

#         group_name = re.match(r'^(Group_\d+)_', file_name).group(1)
#         run_name = re.search(r'_(RUN\d+)_', file_name).group(1)
#         # group_name = re.search(r'_(Group_\d+)_', file_name).group(1)
#         # run_name = re.match(r'^(RUN\d+)_', file_name).group(1)
#         print(group_name,run_name)
#         if group_name == 'Group_5' or group_name == 'Group_2' or group_name == 'Group_4':
#             continue
#         if os.path.isfile(file_path):
#             cur_df=pd.read_csv(file_path,header=0).dropna()
#             cur_df['group']=group_name
#             cur_df['run']=run_name
#             sign_window_df=pd.concat([sign_window_df,cur_df],ignore_index=True)

# sign_window_df['fd'] = sign_window_df['fd'].apply(lambda x: 0 if x < 0 else x)
# sign_window_df = sign_window_df[sign_window_df['fd']<=1]
# sign_window_df=sign_window_df.groupby(['group','start','scaffold']).agg({'fd':'mean','end':'mean'}).reset_index()
# sign_window_df['sign'] = 'Sign'
# # sign_window_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Significant_window_fd_mean.csv',index=False)

if False: #make average fd
    fd_df=pd.DataFrame()
    target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\fd_0504batch'
    for file_name in os.listdir(target_directory):
            if file_name.endswith('SubC_50000_100.csv'):
                fileC_path = os.path.join(target_directory, file_name)
                fileDE_path = os.path.join(target_directory, file_name.replace('SubC', 'SubDE'))
                for i in [fileC_path, fileDE_path]:
                    temp_df = read_fd_df(i)
                    fd_df=fd_df._append(temp_df,ignore_index=True)
    # print(fd_df)
    fd_df['fd'] = fd_df['fd'].apply(lambda x: 0 if x < 0 else x)
    fd_df['fd'] = fd_df['fd'].apply(lambda x: 0 if x > 1 else x)
    fd_df['fd'] = fd_df.apply(lambda x: 0 if x['D'] < 0 else x['fd'], axis=1)
    fd_df=fd_df.groupby(['group','start','scaffold']).agg({'fd':'mean','end':'mean'}).reset_index()
    fd_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Introgressed_50k_window_fd_mean.csv',index=False)

if True: #make average pi
    fd_df=pd.DataFrame()
    target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fst_0504batch'
    for allo_file in os.listdir(target_directory):
            if allo_file.endswith('SubC_50000_100_allo.csv'):
                sym_file = allo_file.replace('allo', 'sym')
                sym_path = os.path.join(target_directory, sym_file)
                sym_df = read_fd_df(sym_path)
                allo_path = os.path.join(target_directory, allo_file)
                allo_df = read_fd_df(allo_path)
                cur_df = pd.concat([allo_df, sym_df[['pi_sympatric','dxy_sympatric_C4', 'Fst_sympatric_C4']]], axis=1)
                # for i in [allo_file, sym_file]:
                #     fileC_path = os.path.join(target_directory, i)
                #     temp_df = read_fd_df(fileC_path)
                fd_df = fd_df._append(cur_df,ignore_index=True)
    print(fd_df)
    fd_df=fd_df.groupby(['group','start','end','scaffold']).agg({'pi_sympatric':'mean','dxy_sympatric_C4':'mean', 'Fst_sympatric_C4':'mean',
                                                                 'pi_allopatric':'mean','dxy_allopatric_C4':'mean', 'Fst_allopatric_C4':'mean'}).reset_index()
    fd_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Introgressed_50k_window_pi_mean.csv',index=False)

# all_df_df=pd.merge(sign_window_df,fd_df,on=['group','start','end','scaffold'],how='outer')
# # print(all_df_df)

# sym_df=pd.DataFrame()
# target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fst'
# for file_name in os.listdir(target_directory):
#     if file_name.endswith('10000_20_sym.csv'):
#         file_path = os.path.join(target_directory, file_name)
#         group_name = re.match(r'^(Group_\d+)_', file_name).group(1)
#         run_name = re.search(r'_(RUN\d+)_', file_name).group(1)
#         # group_name = re.search(r'_(Group_\d+)_', file_name).group(1)
#         # run_name = re.match(r'^(RUN\d+)_', file_name).group(1)
#         print(group_name,run_name)
#         if group_name == 'Group_5' or group_name == 'Group_2' or group_name == 'Group_4':
#             continue
#         if os.path.isfile(file_path):
#             cur_df=pd.read_csv(file_path,header=0,usecols=[0,1,2,5,7,8]).dropna()
#             cur_df['group']=group_name
#             cur_df['run']=run_name
#             sym_df=pd.concat([sym_df,cur_df],ignore_index=True)
# # print(sym_df)
# sym_df = sym_df.groupby(['group','start','scaffold']).agg({'dxy_sympatric_C4':'mean','pi_sympatric':'mean','Fst_sympatric_C4':'mean','end':'mean'}).reset_index()
# fd_pi_df = pd.merge(fd_df,sym_df,on=['group','start','scaffold'],how='outer')
# fd_pi_df = fd_pi_df.dropna(subset=['fd'])
# bins = pd.IntervalIndex.from_breaks(np.arange(0,1,0.015).round(3), closed='left')
# fd_pi_df['fd_group'] = pd.cut(fd_pi_df['fd'], bins)
# # fd_pi_grouped = fd_pi_df.groupby('fd_group')['pi_sympatric'].mean().reset_index()
# fd_pi_df['bin_labels'] = fd_pi_df['fd_group'].astype('category').cat.codes 
# fd_pi_df=pd.merge(fd_pi_df, Chr_subg_df,on='scaffold',how='left')
# # print(fd_pi_df)


# draw_df = pd.merge(all_df_df,sym_df,on=['group','start','scaffold','end'],how='outer')
# draw_df['sign'] = draw_df['sign'].apply(lambda x: 'Nonsign' if pd.isna(x) else x)

# draw_df=draw_df.dropna(axis='index', how='all',subset=['dxy_sympatric_C4','pi_sympatric','Fst_sympatric_C4'])
# draw_df['Fst_sympatric_C4'] = draw_df['Fst_sympatric_C4'].apply(lambda x: 0 if x < 0 else x)
# # print(draw_df)
# fig, axs = plt.subplots(3, 2, figsize=(12, 12))
# sns.boxplot(x='group', y='Fst_sympatric_C4', hue='sign', data=draw_df, palette=['#CA0020', '#0571B0'],
#             fill=False,
#             ax=axs[0,0])
# axs[0,0].set_title('')
# axs[0,0].set_xlabel('Group')
# axs[0,0].set_ylabel('Fst')
# axs[0,0].set_ylim(0,1)
# axs[0,0].grid(True, linestyle='--')

# sns.boxplot(x='group', y='dxy_sympatric_C4', hue='sign', data=draw_df, palette=['#CA0020', '#0571B0'],
#             fill=False,
#             ax=axs[1,0])
# axs[1,0].set_title('')
# axs[1,0].set_xlabel('Group')
# axs[1,0].set_ylabel('Dxy')
# axs[1,0].grid(True, linestyle='--')

# sns.boxplot(x='group', y='pi_sympatric', hue='sign', data=draw_df, palette=['#CA0020', '#0571B0'], 
#             fill=False,
#             ax=axs[2,0])
# axs[2,0].set_title('')
# axs[2,0].set_xlabel('')
# axs[2,0].set_ylabel('Pi')
# axs[2,0].grid(True, linestyle='--')

# for idx,i in enumerate(Chr_subg_df['Subg'].unique()):
#     pi_df_draw_df=fd_pi_df[fd_pi_df['Subg']==i]

#     sns.boxplot(x='fd_group', y='pi_sympatric', data=pi_df_draw_df, width=0.8, color='grey', ax=axs[idx,1])
#     sns.regplot(x='bin_labels', y='pi_sympatric', data=pi_df_draw_df, scatter=False, order=5, color='red', ci=95, ax=axs[idx,1])
#     # sns.regplot(x='fd', y='pi_sympatric', data=fd_pi_df, scatter=False, order=2, color='red', ci=95, ax=axs[3])
#     axs[idx,1].set_title('')
#     axs[idx,1].set_xlabel('fd')
#     axs[idx,1].set_xticks([0, 16, 33, 49, 65])
#     axs[idx,1].set_xticklabels([0, 0.25, 0.5, 0.75, 1.0])
#     axs[idx,1].set_ylabel('Pi')
#     axs[idx,1].grid(True, linestyle='--')

# lines, labels = fig.axes[-2].get_legend_handles_labels()
# fig.legend(lines, labels, loc = 'upper right')

# all_group=draw_df['group'].unique()
# pairs=[]
# for g1 in all_group:
#     set=((g1,'Sign'),(g1,'Nonsign'))
#     pairs.append(set)

# # print(pairs,all_df,all_group)
# anno = Annotator(axs[0,0], pairs=pairs, data=draw_df, x='group', y='Fst_sympatric_C4',hue="sign")
# anno.configure(test='t-test_ind', text_format='star',line_height=0.03,line_width=1,comparisons_correction='bonferroni')
# anno.apply_and_annotate()

# anno = Annotator(axs[1,0], pairs=pairs, data=draw_df, x='group', y='dxy_sympatric_C4',hue="sign")
# anno.configure(test='t-test_ind', text_format='star',line_height=0.03,line_width=1,comparisons_correction='bonferroni')
# anno.apply_and_annotate()

# anno = Annotator(axs[2,0], pairs=pairs, data=draw_df, x='group', y='pi_sympatric',hue="sign")
# anno.configure(test='t-test_ind', text_format='star',line_height=0.03,line_width=1,comparisons_correction='bonferroni')
# anno.apply_and_annotate()
# plt.tight_layout()
# # plt.show()
# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fst_dxy_pi_sign_orno.pdf')