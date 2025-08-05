'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-06 17:08:11
LastEditors: Ne0tea
LastEditTime: 2025-07-02 00:59:28
'''
import os
import pronto
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np

# pd.set_option('display.max_rows', None)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

result_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Significant_test'
go_anno = r'E:\Bio_analysis\Database\go-basic.obo'
pfam_anno = r'E:\Bio_analysis\Database\Pfam-A.clans.tsv'

if True:
    go_ont = pronto.Ontology(go_anno)
    All_term_dic = {}
    for file_name in os.listdir(result_dir):
        if file_name.endswith('_go_sum.fisher'):
            file_path = os.path.join(result_dir, file_name)
            with open(file_path, 'r') as annof:
                for i in  annof:
                    line = i.strip().split('\t')
                    term,all_number,term_in_all,target_number,term_in_target = line[:5]
                    term = term.split('(')[0]
                    if term not in All_term_dic:
                        All_term_dic[term] = [int(all_number),int(term_in_all),int(target_number),int(term_in_target)]
                    else:
                        cur_number = All_term_dic[term]
                        All_term_dic[term] = [int(all_number),int(term_in_all),int(cur_number[2])+int(target_number),int(cur_number[3]) + int(term_in_target)]
    for go in list(All_term_dic.keys()):
        All_term_dic[go_ont.get(go).name] = All_term_dic[go]
        All_term_dic.pop(go)

    term_df = pd.DataFrame()
    for go in All_term_dic:
        cur_pair = All_term_dic[go]
        table = np.array([[cur_pair[3], cur_pair[2]],
                          [cur_pair[1], cur_pair[0]]])
        _, p_value = fisher_exact(table, alternative='greater')
        term_df = term_df._append({'GO': go, 'p_value': p_value, 'Term_count': cur_pair[3]}, ignore_index=True)

    term_df = term_df[term_df['p_value'] < 0.05]
    rejected, corrected_p, _, _ = multipletests(term_df['p_value'], method='fdr_bh')
    term_df['Q_value'] = corrected_p
    term_df['logFDR'] = term_df['Q_value'].apply(lambda x: -np.log10(x))
else:
    term_df = pd.read_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\GO_pfam_0504batch\Group_12_RUN3.gene_Agrigo',sep='\t', header=0)
    term_df['logFDR'] = term_df['FDR'].apply(lambda x: -np.log10(x))
    term_df['Fold'] = term_df.apply(lambda x: (x['queryitem']/x['querytotal']) / (x['bgitem']/x['bgtotal']),axis=1)

go_2_show = term_df[term_df['logFDR'] > 4]
print(go_2_show)
go_2_show = go_2_show[go_2_show['term_type'] == 'P']
print(go_2_show)
# go_2_show['Category'] = np.where(
#     go_2_show['Term'].str.contains('stress|response|abiotic', case=False),  # 匹配关键词（忽略大小写）
#     'Response',
#     np.where(
#         go_2_show['Term'].str.contains('germination|development|formation|elongation', case=False),  # 多关键词用 | 分隔
#         'Develop',
#         'Other'
#     )
# )
# go_2_show = go_2_show[go_2_show['Category'] != 'Other']

# # plt.figure(figsize=(5, 20))
# # plt.hlines(go_2_show['Term'], xmin=0, xmax=go_2_show['logFDR'], color='lightgrey', linewidth=2)
# # sns.scatterplot(data=go_2_show, x="logFDR", y="Term", size = 'queryitem', sizes=(50, 200), hue='Fold', palette='viridis')

# go_p=sns.catplot(data=go_2_show,x="logFDR",y="Term",markersize = 4,
#                  row="Category",kind="point",hue='Fold',sharey=False,palette="viridis",
#                  errorbar=None,height=3,aspect=2,orient='h'
# )
# go_p.set_titles("{row_name}")
# for ax in go_p.axes.flat:
#     ax.grid(True,linestyle="--",alpha=0.5,color="gray",axis="x")
# # go_p.set_grid(True, linestyle='--',axis='x')
# # plt.show()
# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fd_gene_enrichment_GO.pdf') 

if False:
    # if for Pfam
    pfam_anno_dic = {} 
    with open(pfam_anno,'r') as pfamf:
        for i in pfamf:
            line=i.strip().split('\t')
            pfam_anno_dic[line[0]] = line[4]

    All_term_dic = {}
    for file_name in os.listdir(result_dir):
        if file_name.endswith('_pfam_sum.fisher'):
            file_path = os.path.join(result_dir, file_name)
            with open(file_path, 'r') as annof:
                for i in  annof:
                    line = i.strip().split('\t')
                    term,all_number,term_in_all,target_number,term_in_target = line[:5]
                    term = term.split('(')[0]
                    if term not in All_term_dic:
                        All_term_dic[term] = [int(all_number),int(term_in_all),int(target_number),int(term_in_target)]
                    else:
                        cur_number = All_term_dic[term]
                        if int(cur_number[3]) > int(term_in_target):
                            # All_term_dic[term] = [int(all_number),int(term_in_all),int(cur_number[2])+int(target_number),int(cur_number[3]) + int(term_in_target)]
                            All_term_dic[term] = [int(all_number),int(term_in_all),int(cur_number[2]),int(cur_number[3])]
                        else:
                            All_term_dic[term] = [int(all_number),int(term_in_all),int(target_number),int(term_in_target)]

    term_df = pd.DataFrame()
    for go in All_term_dic:
        cur_pair = All_term_dic[go]
        table = np.array([[cur_pair[3], cur_pair[2]],
                          [cur_pair[1], cur_pair[0]]])
        _, p_value = fisher_exact(table, alternative='greater')
        term_df = term_df._append({'GO': go, 'p_value': p_value, 'Term_count': cur_pair[3], 'Fold':(cur_pair[3]/cur_pair[2]) / (cur_pair[1]/cur_pair[0])}, ignore_index=True)

    rejected, corrected_p, _, _ = multipletests(term_df['p_value'], method='fdr_bh')
    term_df['Q_value'] = corrected_p
    term_df['Name'] = term_df['GO'].apply(lambda x: pfam_anno_dic[x] if x in pfam_anno_dic else x)
    pfam_2_show = term_df[(term_df['p_value'] < 0.01) & (term_df['Term_count'] > 10) & (~term_df['Name'].str.contains('unknown')) \
                        & (term_df['GO'] != 'PF07732') & (term_df['GO'] != 'PF00394')]

    pfam_2_show = pfam_2_show.sort_values('p_value', ascending=False)
    pfam_2_show['LogP'] = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))
    print(pfam_2_show)

    labels = pfam_2_show['Name']
    values = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))
    Folds = pfam_2_show['Fold']

    plt.figure(figsize=(5, 6.5))
    sns.set_style("white")
    sns.despine(top=True, right=True)
    plt.hlines(labels, xmin=0, xmax=values, color='lightgrey', linewidth=2)
    sns.scatterplot(data=pfam_2_show, x="LogP", y="Name", size = 'Term_count', sizes=(75, 250),
                    hue='Fold', palette='coolwarm', hue_norm=Normalize(vmin=1, vmax=4.5))

    plt.xticks(ticks=[0,5,10])
    plt.yticks(fontsize=6)
    # plt.grid(True, linestyle='--',axis='x')
    plt.ylabel('')
    # plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fd_gene_enrichment_pfam.pdf') 