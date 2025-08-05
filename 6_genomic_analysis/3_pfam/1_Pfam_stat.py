'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-28 21:54:04
LastEditors: Ne0tea
LastEditTime: 2025-05-10 21:30:04
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

pd.set_option('display.max_rows', None)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
diff_class_color = {'Digitaria':'#3a4567', 'Echinochloa':'#98cae8', 'Triticum':'#e39466'}

def read_space_file(file, pre):
    result_df = pd.DataFrame()
    with open(file, 'r') as f:
        for i in f:
            line=i.strip().split(' ',1)
            result_df = result_df._append({pre:int(line[0]), 'Pfam':line[1]}, ignore_index=True)
    return result_df

def get_merge_stat_df(df1, df2, df3):
    df1_pre, df2_pre, df3_pre = os.path.basename(df1).split('_')[0], os.path.basename(df2).split('_')[0], os.path.basename(df3).split('_')[0]
    # drad_df = read_space_file(df1, df1_pre)
    # dmil_df = read_space_file(df2, df2_pre)
    # dsan_df = read_space_file(df3, df3_pre)

    drad_df  = pd.read_csv(df1,sep='\s+', names=[df1_pre, 'Pfam'])
    dmil_df  = pd.read_csv(df2,sep='\s+', names=[df2_pre, 'Pfam'])
    dmil_df = dmil_df.groupby('Pfam').sum()
    dsan_df  = pd.read_csv(df3,sep='\s+', names=[df3_pre, 'Pfam'])
    dsan_df = dsan_df.groupby('Pfam').sum()

    all_tf_df = pd.merge(drad_df, dmil_df, on='Pfam', how='outer')
    all_tf_df = pd.merge(all_tf_df, dsan_df, on='Pfam', how='outer')
    all_tf_df.set_index('Pfam', inplace=True)

    all_tf_df['E'] = all_tf_df[df3_pre] - all_tf_df[df2_pre] - all_tf_df[df1_pre]
    all_tf_df[df1_pre] = all_tf_df[df1_pre] / all_tf_df[df1_pre].sum()
    all_tf_df[df2_pre] = all_tf_df[df2_pre] / all_tf_df[df2_pre].sum()
    all_tf_df[df3_pre] = all_tf_df[df3_pre] / all_tf_df[df3_pre].sum()
    all_tf_df['A'] = (all_tf_df[df1_pre] + all_tf_df[df2_pre]) / (all_tf_df[df1_pre].sum() + all_tf_df[df2_pre].sum())
    all_tf_df['E_r'] = all_tf_df[df3_pre] - all_tf_df['A'] 
    return all_tf_df


def draw_strip_plot(plot_df, anno_dic):
    plot_trans_df = pd.DataFrame()
    D_df = plot_df.loc[:, plot_df.columns.str.endswith('_D')].reset_index()
    D_df['Class'] = 'Digitaria'
    D_df.columns = ['Pfam','Expansion', 'AncestorN', 'Expansion_ratio', 'Class']
    E_df = plot_df.loc[:, plot_df.columns.str.endswith('_E')].reset_index()
    E_df['Class'] = 'Echinochloa'
    E_df.columns = ['Pfam','Expansion', 'AncestorN', 'Expansion_ratio', 'Class']
    T_df = plot_df.loc[:, plot_df.columns.str.endswith('_T')].reset_index()
    T_df['Class'] = 'Triticum'
    T_df.columns = ['Pfam','Expansion', 'AncestorN', 'Expansion_ratio', 'Class']
    for i in D_df, E_df, T_df:
        plot_trans_df = pd.concat([plot_trans_df, i], axis=0).reset_index(drop=True)
    plot_trans_df['Name'] = plot_trans_df['Pfam'].apply(lambda x: anno_dic[x] if x in anno_dic else x)
    plt.figure(figsize=(8, 4))
    sns.stripplot(x="Expansion_ratio",y="Name", hue = 'Class', data=plot_trans_df, 
                  palette=diff_class_color, jitter=False, orient='h')
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    plt.show()


pfam_anno = r'E:\Bio_analysis\Database\Pfam-A.clans.tsv'
pfam_anno_dic = {}
with open(pfam_anno,'r') as pfamf:
    for i in pfamf:
        line=i.strip().split('\t')
        pfam_anno_dic[line[0]] = line[4]

Drad_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Drad_pfam.summary'
Dmil_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Dmil_pfam.summary'
Dsan_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Dsan_pfam.summary'
Digitaria_pfam_df = get_merge_stat_df(Drad_tf_summary_file, Dmil_tf_summary_file, Dsan_tf_summary_file)

Ehap_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Ehap_pfam.summary'
Eory_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Eory_pfam.summary'
Ecru_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Ecru_pfam.summary'
Echinochloa_pfam_df = get_merge_stat_df(Ehap_tf_summary_file, Eory_tf_summary_file, Ecru_tf_summary_file)

Atau_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Atau_pfam.summary'
Ttur_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Ttur_pfam.summary'
Taes_tf_summary_file= r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Taes_pfam.summary'
Triticum_pfam_df = get_merge_stat_df(Atau_tf_summary_file, Ttur_tf_summary_file, Taes_tf_summary_file)

Dsan_contract_list = Digitaria_pfam_df[Digitaria_pfam_df['E_r']<=-0.0002].index.to_list()
Dsan_expansion_list = Digitaria_pfam_df[Digitaria_pfam_df['E_r']>=0.0002].index.to_list()

All_df = pd.merge(Digitaria_pfam_df, Echinochloa_pfam_df, left_index=True, right_index=True, how='left',suffixes=('_D','_E'))
All_df = pd.merge(All_df, Triticum_pfam_df,left_index=True, right_index=True, how='left')
All_df = All_df.rename(columns={'E':'E_T', 'A':'A_T', 'E_r':'E_r_T'})
All_df['Name'] = All_df.index.map(lambda x: pfam_anno_dic[x] if x in pfam_anno_dic else x)

All_contract_df = All_df.loc[Dsan_contract_list, All_df.columns.str.endswith(('_D', '_E', '_T', 'Name'))]
All_contract_df = All_contract_df[~All_contract_df['Name'].str.contains('unknown')]
All_contract_df = All_contract_df.loc[All_contract_df.sort_values('E_r_D', ascending=False).index.tolist()]
All_expansion_df = All_df.loc[Dsan_expansion_list, All_df.columns.str.endswith(('_D', '_E', '_T', 'Name'))]
All_expansion_df = All_expansion_df[~All_expansion_df['Name'].str.contains('unknown')]
All_expansion_df = All_expansion_df.loc[All_expansion_df.sort_values('E_r_D', ascending=False).index.tolist()]
# print(All_expansion_df)
draw_strip_plot(All_expansion_df, pfam_anno_dic)
draw_strip_plot(All_contract_df, pfam_anno_dic)
# df.to_parquet("data.parquet", engine='pyarrow')
# # 读取
# df = pd.read_parquet("data.parquet")