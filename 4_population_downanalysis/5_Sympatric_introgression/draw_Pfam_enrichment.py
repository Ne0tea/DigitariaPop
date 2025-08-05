'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-06 17:08:11
LastEditors: Ne0tea
LastEditTime: 2025-07-07 19:42:25
'''
import os
import pronto
import pandas as pd
import seaborn as sns
import PyComplexHeatmap as pch
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np

pd.set_option('display.max_rows', None)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

result_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Significant_test'
go_anno = r'E:\Bio_analysis\Database\go-basic.obo'
pfam_anno = r'E:\Bio_analysis\Database\Pfam-A.clans.tsv'
pop_name_convert_dic={'Group_12':'JN','Group_1':'LF','Group_2':'HZ','Group_3':'SX','Group_4':'DZ','Group_6':'SJZ','Group_7':'ZJ'}

def run_fisher_and_filter(term_dic, term_count=10):

    term_df = pd.DataFrame()
    for go in term_dic:
        cur_pair = term_dic[go]
        table = np.array([[cur_pair[3], cur_pair[2]],
                          [cur_pair[1], cur_pair[0]]])
        _, p_value = fisher_exact(table, alternative='greater')
        term_df = term_df._append({'GO': go, 'p_value': p_value, 'Term_count': cur_pair[3], 'Ref_count': cur_pair[1], 'Fold':(cur_pair[3]/cur_pair[2]) / (cur_pair[1]/cur_pair[0])}, ignore_index=True)

    term_df = term_df[term_df['p_value'] < 0.05]
    rejected, corrected_p, _, _ = multipletests(term_df['p_value'], method='fdr_bh')
    term_df['Q_value'] = corrected_p
    term_df['Name'] = term_df['GO'].apply(lambda x: pfam_anno_dic[x] if x in pfam_anno_dic else x)
    pfam_2_show = term_df[(term_df['p_value'] < 0.05) & (term_df['Term_count'] > term_count) & (~term_df['Name'].str.contains('unknown')) \
                        & (term_df['GO'] != 'PF07732') & (term_df['GO'] != 'PF00394')]
    return pfam_2_show

def draw_overall_enrichment(pfam_2_show):
    pfam_2_show = pfam_2_show.sort_values('Q_value', ascending=True).head(12)
    print(pfam_2_show)
    labels = pfam_2_show['Name']
    values = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))
    Folds = pfam_2_show['Fold']
    pfam_2_show['LogP'] = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))
    plt.figure(figsize=(3, 5))
    sns.set_style("white")
    sns.despine(top=True, right=True)
    sns.scatterplot(data=pfam_2_show, x="LogP", y="Name", size = 'Term_count', sizes=(75, 250),
                    hue='Fold', palette='vlag', hue_norm=Normalize(vmin=2, vmax=12))
    plt.hlines(labels, xmin=0, xmax=values, color='lightgrey', linewidth=3)

    plt.xticks(ticks=[0,5,10])
    plt.yticks(fontsize=12)
    # plt.grid(True, linestyle='--',axis='x')
    plt.ylabel('')
    # plt.tight_layout()
    plt.legend(loc="upper right")
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fd_gene_enrichment_pfam.pdf') 

def draw_per_enrichment(pfam_2_show, file_prefix):
    labels = pfam_2_show['Name']
    values = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))
    Folds = pfam_2_show['Fold']
    pfam_2_show = pfam_2_show.sort_values('p_value', ascending=False)
    pfam_2_show['LogP'] = pfam_2_show['p_value'].apply(lambda x: -np.log10(x))

    plt.figure(figsize=(5, 4))
    sns.set_style("white")
    sns.despine(top=True, right=True)
    sns.scatterplot(data=pfam_2_show, x="LogP", y="Name", size = 'Term_count', sizes=(75, 250),
                    hue='Fold', palette='coolwarm', hue_norm=Normalize(vmin=1, vmax=4.5))
    plt.hlines(labels, xmin=0, xmax=values, color='lightgrey', linewidth=2)

    plt.xticks(ticks=[0,5,10])
    plt.yticks(fontsize=12)
    # plt.grid(True, linestyle='--',axis='x')
    plt.ylabel('')
    # plt.tight_layout()
    plt.legend(loc="upper right")
    # plt.show()
    plt.savefig(rf'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Fd_GE_pfam_{file_prefix}.pdf') 

def draw_dotcomplex(data1, target_pfam_list, pop_list):
    # rename_PF_dic = {'NB-ARC domain':'NB-ARC', 'D-mannose binding lectin':'D-mannose binding lectin', 'No apical meristem-associated C-terminal domain':'NAM C-terminal',
    #                  'UDP-glucoronosyl and UDP-glucosyl transferase':'UDP-glucoronosyl and UDP-glucosyl transferase',
    #                  'AP2 domain':'AP2', 'GRAS domain family':'GRAS', 
    #                  'ABC transporter':'ABC transporter','Cytochrome P450':'P450', 'Hsp70 protein':'Hsp70',
    #                  'FAR1 DNA-binding domain':'FAR1'}
    # PF_order = ['NB-ARC', 'D-mannose binding lectin', 'UDP-glucoronosyl and UDP-glucosyl transferase', 'NAM C-terminal',
    #             'AP2', 'GRAS', 'ABC transporter', 'P450', 'Hsp70','FAR1']
    # data1['Pfam'] = data1['Pfam'].map(rename_PF_dic)

    # col_ha = pch.HeatmapAnnotation(
    #                            Subfamily=pch.anno_simple(sp_df.Subfamily,
    #                                                      legend=False,add_text=True,
    #                                                      height=3,
    #                                                      colors=subfamily_class_colos),
    #                            Subtribe=pch.anno_simple(sp_df.Subtribe,
    #                                                     legend=False,add_text=True,
    #                                                     height=3,
    #                                                     colors=subtribe_class_colos),
    #                            verbose=0,label_side='right',label_kws={'horizontalalignment':'left'})
    sp_df = pd.DataFrame(pop_list)
    sp_df.columns = ['Pop']
    sp_df['dd'] = sp_df['Pop']
    sp_df.set_index('dd', inplace=True)
    col_ha = pch.HeatmapAnnotation(
                               Subfamily=pch.anno_simple(sp_df.Pop,
                                                         legend=False,add_text=True,
                                                         height=3),
                               verbose=0,label_side='right',label_kws={'horizontalalignment':'left'})

    plt.figure(figsize=(4, 6))

    cm = pch.DotClustermapPlotter(data=data1, x='Pop',y='Name',value='logQ', c='logQ',s='Fold',
                                #   zscore=True,
                                  top_annotation = col_ha,
                                  show_rownames=True,show_colnames=True,
                                  col_names_side='Top',row_names_side='left',
                                  row_cluster=False,
                                  col_cluster=False,
                                #   xticklabels=True,
                                  cmap='coolwarm',
                                  vmax=2.4,
                                  vmin=1.3,
                                #   col_split=sp_df.Subtribe,
                                #   col_split_gap = 1
                                  x_order=pop_list,
                                  y_order=target_pfam_list
                             )
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Introgressed_pop_pfam.pdf', format='pdf')

def draw_upset(long_data):
    from upsetplot import from_memberships
    from upsetplot import UpSet

    set_list = long_data.groupby("Name")["Pop"].apply(tuple).to_list()

    set_counts_series = long_data.groupby("Name")["Pop"].apply(tuple).value_counts()
    tuple_count = []
    for i in set_list:
        tuple_count.append(set_counts_series[i])

    upset_data = from_memberships(set_list, data=tuple_count)
    value_count_df = pd.DataFrame(long_data.groupby("Name")["Pop"].apply(tuple).value_counts())
    shared_pfam_count_over2 = 0
    shared_pfam_count_over3 = 0
    shared_pfam_count_over4 = 0
    for set_pop, count_df in value_count_df.iterrows():
        cur_set_compo = len(list(set_pop))
        if cur_set_compo == 2 :
            shared_pfam_count_over2 += count_df['count']
        elif cur_set_compo == 3 :
            shared_pfam_count_over3 += count_df['count']
        elif cur_set_compo >= 4 :
            shared_pfam_count_over4 += count_df['count']
    print(shared_pfam_count_over2)
    print(shared_pfam_count_over3)
    print(shared_pfam_count_over4)

    # upset = UpSet(upset_data, subset_size="count", show_counts=True, element_size=20)

    # upset.plot()
    # plt.title("")
    # # plt.show()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Introgressed_window\Introgressed_pop_pfam_upset.pdf', format='pdf')

if True:
    target_pfam_list = ['ABC transporter','Myb-like DNA-binding domain', 'UDP-glucoronosyl and UDP-glucosyl transferase',
                        'FBD', 'FAR1 DNA-binding domain',
                        'Cytochrome P450', 'C2 domain', 'AP2 domain',
                        'Hsp20/alpha crystallin family',
                        'D-mannose binding lectin']
    pop_list = [ 'DZ', 'SJZ', 'JN', 'HZ', 'ZJ', 'LF', 'SX']
    # if for Pfam
    pfam_anno_dic = {}
    with open(pfam_anno,'r') as pfamf:
        for i in pfamf:
            line=i.strip().split('\t')
            pfam_anno_dic[line[0]] = line[4]

    All_term_dic = {}
    single_pop_term_dic = {}
    for file_name in os.listdir(result_dir):
        if file_name.endswith('_pfam_sum.fisher'):
            popname = pop_name_convert_dic[file_name.rsplit('_',6)[0]]
            file_path = os.path.join(result_dir, file_name)
            single_pop_term_dic[popname] = {}
            with open(file_path, 'r') as annof:
                for i in  annof:
                    line = i.strip().split('\t')
                    term,all_number,term_in_all,target_number,term_in_target = line[:5]
                    term = term.split('(')[0]
                    single_pop_term_dic[popname][term] = [int(all_number),int(term_in_all),int(target_number),int(term_in_target)]
                    if term not in All_term_dic:
                        All_term_dic[term] = [int(all_number),int(term_in_all),int(target_number),int(term_in_target)]
                    else:
                        cur_number = All_term_dic[term]
                        All_term_dic[term] = [int(all_number),int(term_in_all),
                                              int(cur_number[2]) + int(target_number), 
                                              int(cur_number[3]) + int(term_in_target)]
    All_pfam_2_show = pd.DataFrame()
    for i in single_pop_term_dic:
        pfam_2_show = run_fisher_and_filter(single_pop_term_dic[i], 3)
        pfam_2_show['Pop'] = i
        All_pfam_2_show = pd.concat([All_pfam_2_show, pfam_2_show])
    All_pfam_2_show.reset_index(drop=True, inplace=True)
    All_pfam_2_show = All_pfam_2_show[All_pfam_2_show['GO']!='PF13417']
    All_pfam_2_show.to_csv('temp.csv',index=False)
    # All_pfam_2_show_wide = All_pfam_2_show.pivot(
    #     index='Name',
    #     columns='Pop',
    #     values='Q_value'
    # )
    # print(All_pfam_2_show_wide)
    # draw_upset(All_pfam_2_show)

    #draw PFAM dot heatmap
    All_pfam_2_show = All_pfam_2_show[All_pfam_2_show['Name'].isin(target_pfam_list)]
    All_pfam_2_show['logQ'] = All_pfam_2_show['Q_value'].apply(lambda x: -np.log10(x))
    print(All_pfam_2_show)
    # draw_dotcomplex(All_pfam_2_show, target_pfam_list, pop_list)

    # xxx = All_pfam_2_show.groupby('Name').sum()
    # print(xxx)
    # All_pfam_2_show.to_csv('temp.csv',index=False)
    # pfam_2_show = run_fisher_and_filter(All_term_dic)
    # draw_overall_enrichment(pfam_2_show)
