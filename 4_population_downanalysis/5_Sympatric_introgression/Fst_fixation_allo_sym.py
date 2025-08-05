'''
Descripttion: 
Author: Ne0tea
version: env path: 
Date: 2024-09-20 22:39:30
LastEditors: Ne0tea
LastEditTime: 2025-07-25 22:55:22
'''
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from statannotations.Annotator import Annotator

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
pd.set_option('display.max_rows', None)

genome_base=1347733405
allo_df=pd.DataFrame()
sym_df=pd.DataFrame()
target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fst_0504batch'
for file_name in os.listdir(target_directory):
    if file_name.endswith('SubC_100000_200_allo.csv'):
        file_path = os.path.join(target_directory, file_name)
        group_name = re.search(r'_(Group_\d+)_', file_name).group(1)
        run_name = re.match(r'^(RUN\d+)_', file_name).group(1)
        if os.path.isfile(file_path):
            cur_df=pd.read_csv(file_path,header=0,usecols=[0,1,2,3,8])
            cur_df = cur_df.fillna(0)
            cur_df['group']=group_name
            cur_df['run']=run_name
            cur_df = cur_df.rename(columns={'Fst_allopatric_C4': 'Fst'})
            # cur_df = cur_df[cur_df['Fst'] >= 0]
            # cur_df['Fst'] = cur_df['Fst'].apply(lambda x: 0 if x < 0 else x)
            allo_df=pd.concat([allo_df,cur_df],ignore_index=True)

    elif file_name.endswith('SubC_100000_200_sym.csv'):
        file_path = os.path.join(target_directory, file_name)

        group_name = re.search(r'_(Group_\d+)_', file_name).group(1)
        run_name = re.match(r'^(RUN\d+)_', file_name).group(1)
        if os.path.isfile(file_path):
            cur_df=pd.read_csv(file_path,header=0,usecols=[0,1,2,3,8])
            cur_df = cur_df.fillna(0)
            cur_df['group']=group_name
            cur_df['run']=run_name
            cur_df = cur_df.rename(columns={'Fst_sympatric_C4': 'Fst'})
            # cur_df = cur_df[cur_df['Fst'] >= 0]
            # cur_df['Fst'] = cur_df['Fst'].apply(lambda x: 0 if x < 0 else x)
            sym_df=pd.concat([sym_df,cur_df],ignore_index=True)

# allo_df = allo_df.groupby(['group', 'scaffold', 'mid']).agg({'Fst':'mean'}).reset_index()
# sym_df = sym_df.groupby(['group', 'scaffold', 'mid']).agg({'Fst':'mean'}).reset_index()
allo_df['location'] = 'Allopatric'
sym_df['location'] = 'Sympatric'
all_df = pd.merge(allo_df, sym_df, on=['group', 'start', 'end', 'scaffold', 'mid', 'run'], how='inner')
allo_inner_df = all_df[['scaffold','start','end','mid','Fst_x','group','run','location_x']].rename(columns={'Fst_x': 'Fst','location_x':'location'})
sym_inner_df = all_df[['scaffold','start','end','mid','Fst_y','group','run','location_y']].rename(columns={'Fst_y': 'Fst','location_y':'location'})
all_df = pd.concat([allo_inner_df, sym_inner_df], ignore_index=True)
all_group=all_df['group'].unique()
print(all_df[all_df['group'] == 'Group_2'].groupby('location')['Fst'].mean())

plt.figure(figsize=(8, 3))
ax = sns.boxplot(x='group', y='Fst', hue='location', data=all_df,
                 showfliers=False,
                 fill=False,
                 width=0.6,
                 linewidth=1.5,
                 palette=['#92C5DE', '#0571B0'])
plt.ylim(0,1)
plt.xlabel('')
plt.grid(True, linestyle='--')
pairs=[]
for g1 in all_group:
    set=((g1,'Allopatric'),(g1,'Sympatric'))
    pairs.append(set)

# print(pairs,all_df,all_group)
anno = Annotator(ax, pairs=pairs, data=all_df, x='group', y='Fst',hue="location")
anno.configure(test='t-test_ind', text_format='star',line_height=0.03,line_width=1,comparisons_correction='bonferroni')
anno.apply_and_annotate()

plt.legend()
plt.show()
# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Fst_distribution_in_allo_sym.pdf')