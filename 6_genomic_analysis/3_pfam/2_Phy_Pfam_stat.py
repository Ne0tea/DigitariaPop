'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-10 23:22:42
LastEditors: Ne0tea
LastEditTime: 2025-07-31 15:45:28
'''
import os
import pandas as pd
import numpy as np
import PyComplexHeatmap as pch
import matplotlib.pylab as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
target_PF = ['PF03101','PF00931', 'PF01453', 'PF14303' ,'PF00847', 'PF00012', 'PF00201', 'PF00005', 'PF00067', 'PF03514', 'PF05686', 'PF00248']
subtribe_class = {'Osat':'Out','Eind':'Out','Ehap':'Boivinellinae','Ecru-A':'Boivinellinae', 'Ecru-B':'Boivinellinae',
                  'Ecru-C':'Boivinellinae', 'Eory-A':'Boivinellinae', 'Eory-B':'Boivinellinae', 'Sita':'Cenchrinae',
                  'Phal':'Panicinae','Dexi-B':'Anthephorinae','Dexi-A':'Anthephorinae','Dsan-C':'Anthephorinae',
                  'Drad':'Anthephorinae','Dsan-E':'Anthephorinae','Dmil-E':'Anthephorinae','Dsan-D':'Anthephorinae','Dmil-D':'Anthephorinae'}
subtribe_class_colos = {'Out':'#707070', 'Boivinellinae':'#606c38', 'Cenchrinae':'#ccd5ae', 'Panicinae':'#faedcd', 'Anthephorinae':'#dda15e'}
subfamily_class = {'Osat':'ORY','Eind':'CHL','Ehap':'PAN','Ecru-A':'PAN', 'Ecru-B':'PAN',
                  'Ecru-C':'PAN', 'Eory-A':'PAN', 'Eory-B':'PAN', 'Sita':'PAN',
                  'Phal':'PAN','Dexi-B':'PAN','Dexi-A':'PAN','Dsan-C':'PAN',
                  'Drad':'PAN','Dsan-E':'PAN','Dmil-E':'PAN','Dsan-D':'PAN','Dmil-D':'PAN'}
subfamily_class_colos = {"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1"}

colors = ["#006d77", "#e5eff2", "#e29578"]
custom_cmap = LinearSegmentedColormap.from_list("Costum", colors)
plt.register_cmap(name="Costum", cmap=custom_cmap)

def read_phy_order(config:str):
    clade_subg_dic={}

    with open(config,'r') as c_file:
        phlo=0
        subg=0
        sp=0
        for i in c_file:
            if i.startswith('>'):
                if i.startswith('>phylogeny'):
                    phlo,subg=1,0
                if i.startswith('>subg list'):
                    phlo,subg=0,1
            else:
                if phlo:
                    sub_phlo=i.strip()
                if subg:
                    line=i.strip().split()
                    clade_subg_dic[line[0]]=line[1:]
    return clade_subg_dic

def draw_dotcomplex(data1):
    species_order = ['Osat','Eind','Ehap','Ecru-A', 'Ecru-B','Ecru-C', 'Eory-A', 'Eory-B', 'Sita','Phal','Dexi-B','Dexi-A','Dsan-C','Drad','Dsan-E','Dmil-E','Dsan-D','Dmil-D']
    sp_df = pd.DataFrame(species_order)
    sp_df.columns = ['Species']

    rename_PF_dic = {'NB-ARC domain':'NB-ARC', 'D-mannose binding lectin':'D-mannose binding lectin', 'No apical meristem-associated C-terminal domain':'NAM C-terminal',
                     'UDP-glucoronosyl and UDP-glucosyl transferase':'UDP-glucoronosyl and UDP-glucosyl transferase',
                     'AP2 domain':'AP2', 'GRAS domain family':'GRAS', 
                     'ABC transporter':'ABC transporter','Cytochrome P450':'P450', 'Hsp70 protein':'Hsp70',
                     'FAR1 DNA-binding domain':'FAR1'}
    PF_order = ['NB-ARC', 'D-mannose binding lectin', 'UDP-glucoronosyl and UDP-glucosyl transferase', 'NAM C-terminal',
                'AP2', 'GRAS', 'ABC transporter', 'P450', 'Hsp70','FAR1']
    data1['Pfam'] = data1['Pfam'].map(rename_PF_dic)
    sp_df['Subtribe'] = sp_df['Species'].map(subtribe_class)
    sp_df['Subfamily'] = sp_df['Species'].map(subfamily_class)
    sp_df = sp_df.set_index('Species')
    col_ha = pch.HeatmapAnnotation(
                               Subfamily=pch.anno_simple(sp_df.Subfamily,
                                                         legend=False,add_text=True,
                                                         height=3,
                                                         colors=subfamily_class_colos),
                               Subtribe=pch.anno_simple(sp_df.Subtribe,
                                                        legend=False,add_text=True,
                                                        height=3,
                                                        colors=subtribe_class_colos),
                               verbose=0,label_side='right',label_kws={'horizontalalignment':'left'})

    plt.figure(figsize=(6, 4))
    # print(data1)
    cm = pch.DotClustermapPlotter(data=data1, x='Species',y='Pfam',value='zscore', c='zscore',s=0.5,
                                #   zscore=True,
                                  cmap='Costum',
                                  top_annotation = col_ha,
                                  show_rownames=True,show_colnames=True,
                                  col_names_side='bottom',row_names_side='left',
                                  row_cluster=False,col_cluster=False,
                                  vmax=2.5,vmin=-2.5,
                                  col_split=sp_df.Subtribe,
                                  col_split_order = list(subtribe_class_colos.keys()),
                                  col_split_gap = 1,
                                  x_order=species_order,
                                  y_order=PF_order
                             )
    plt.show()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Pfam_stat_modi.pdf', format='pdf')

pfam_anno = r'E:\Bio_analysis\Database\Pfam-A.clans.tsv'
pfam_anno_dic = {}
with open(pfam_anno,'r') as pfamf:
    for i in pfamf:
        line=i.strip().split('\t')
        pfam_anno_dic[line[0]] = line[4]

phylogeny_class = read_phy_order(r'E:\Bio_analysis\HGT_newcollection\2_HGTfinder_workpath\Poaceae_config.txt')
sp_ord = phylogeny_class['PAN'] + ['Drad', 'Dmil-D', 'Dmil-E']
target_directory = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\General_pfam'
phy_pfam_df = pd.DataFrame()
all_pfam_number_dic = {}
for file_name in os.listdir(target_directory):
    cur_pfam_file = os.path.join(target_directory,file_name)
    cur_sp_name = file_name.split('_')[0]
    sp_pfam_df  = pd.read_csv(cur_pfam_file,sep='\s+', names=[cur_sp_name, 'Pfam'])
    all_pfam_number_dic[cur_sp_name] = sp_pfam_df[cur_sp_name].sum()
    if phy_pfam_df.empty:
        phy_pfam_df = sp_pfam_df
    else:
        phy_pfam_df = pd.merge(phy_pfam_df, sp_pfam_df, on='Pfam', how='outer')
phy_pfam_df = phy_pfam_df.set_index('Pfam')

all_sp_list = phy_pfam_df.columns.to_list()
max_cols = phy_pfam_df.idxmax(axis=1)

# result = phy_pfam_df[max_cols.str.contains('Dsan')].reset_index()
# result['Name'] = result['Pfam'].apply(lambda x: pfam_anno_dic[x] if x in pfam_anno_dic else x)
# result.to_clipboard(index=False)

phy_pfam_df['Name'] =  phy_pfam_df.index.map(pfam_anno_dic)
ord_pfam_df = phy_pfam_df.loc[phy_pfam_df[sp_ord].min(axis=1) > 5, sp_ord].apply(stats.zscore, axis=1)
ord_pfam_df['Name'] =  ord_pfam_df.index.map(pfam_anno_dic)
# ord_pfam_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Pfam_analysis\Pfam_compare_stat.csv')
plot_pfam_df = phy_pfam_df.loc[target_PF].reset_index()
plot_pfam_df['Pfam'] = plot_pfam_df['Pfam'].apply(lambda x: pfam_anno_dic[x] if x in pfam_anno_dic else x)

plot_pfam_long_df = pd.melt(
    plot_pfam_df,
    id_vars=["Pfam"],
    value_vars=all_sp_list,
    var_name="Species",
    value_name="count"
)
plot_pfam_long_df['zscore'] = plot_pfam_long_df.groupby('Pfam')['count'].transform(
    lambda x: (x - x.mean()) / x.std()
)
plot_pfam_long_df['relative_count'] = plot_pfam_long_df.apply(lambda x: x['count'] / all_pfam_number_dic[x['Species']], axis=1)
# print(plot_pfam_long_df)
# draw_dotcomplex(plot_pfam_long_df)


# for domain in ['D-mannose binding lectin', 'NB-ARC domain']:
#     print(domain)
#     D_target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].str.startswith('D'))]['count'].tolist()
#     D_medians = np.median(D_target_pfam_values)
#     D_q1 = np.percentile(D_target_pfam_values, 10)
#     D_q3 = np.percentile(D_target_pfam_values, 90)
#     D_stds = np.std(D_target_pfam_values)
#     print(D_medians, D_q1, D_q3, D_stds)

#     target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].isin(['Phal','Osat','Sita']))]['count'].tolist()
#     medians = np.median(target_pfam_values)
#     q1 = np.percentile(target_pfam_values, 10)
#     q3 = np.percentile(target_pfam_values, 90)
#     stds = np.std(target_pfam_values)
#     print(medians, q1, q3, stds)

#     t_stat, p_value = stats.ttest_ind(D_target_pfam_values, target_pfam_values)
#     print(t_stat, p_value)

# for domain in ['UDP-glucoronosyl and UDP-glucosyl transferase']:
#     print(domain)
#     D_target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].str.startswith('D'))]['count'].tolist()
#     D_medians = np.median(D_target_pfam_values)
#     D_q1 = np.percentile(D_target_pfam_values, 10)
#     D_q3 = np.percentile(D_target_pfam_values, 90)
#     D_stds = np.std(D_target_pfam_values)
#     print(D_medians, D_q1, D_q3, D_stds)

#     target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].isin(['Ehap','Eory-A','Eory-B','Ecru-A','Ecru-B','Ecru-C']))]['count'].tolist()
#     medians = np.median(target_pfam_values)
#     q1 = np.percentile(target_pfam_values, 10)
#     q3 = np.percentile(target_pfam_values, 90)
#     stds = np.std(target_pfam_values)
#     print(medians, q1, q3, stds)

#     t_stat, p_value = stats.ttest_ind(D_target_pfam_values, target_pfam_values)
#     print(t_stat, p_value)

# for domain in ['GRAS domain family', 'UDP-glucoronosyl and UDP-glucosyl transferase', 'AP2 domain']:
#     print(domain)
#     D_target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].str.startswith('Dsan'))]['count'].tolist()
#     D_medians = np.median(D_target_pfam_values)
#     D_q1 = np.percentile(D_target_pfam_values, 10)
#     D_q3 = np.percentile(D_target_pfam_values, 90)
#     D_stds = np.std(D_target_pfam_values)
#     print(D_medians, D_q1, D_q3, D_stds)

#     target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].isin(['Ehap','Eory-A','Eory-B','Ecru-A','Ecru-B','Ecru-C']))]['count'].tolist()
#     medians = np.median(target_pfam_values)
#     q1 = np.percentile(target_pfam_values, 10)
#     q3 = np.percentile(target_pfam_values, 90)
#     stds = np.std(target_pfam_values)
#     print(medians, q1, q3, stds)

#     t_stat, p_value = stats.ttest_ind(D_target_pfam_values, target_pfam_values)
#     print(t_stat, p_value)

# for domain in ['FAR1 DNA-binding domain']:
#     print(domain)
#     D_target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].str.startswith('Dsan'))]['count'].tolist()
#     D_medians = np.median(D_target_pfam_values)
#     D_q1 = np.percentile(D_target_pfam_values, 10)
#     D_q3 = np.percentile(D_target_pfam_values, 90)
#     D_stds = np.std(D_target_pfam_values)
#     print(D_medians, D_q1, D_q3, D_stds)

#     target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].isin(['Dmil-D','Dmil-E','Drad']))]['count'].tolist()
#     medians = np.median(target_pfam_values)
#     q1 = np.percentile(target_pfam_values, 10)
#     q3 = np.percentile(target_pfam_values, 90)
#     stds = np.std(target_pfam_values)
#     print(medians, q1, q3, stds)

#     t_stat, p_value = stats.ttest_ind(D_target_pfam_values, target_pfam_values)
#     print(t_stat, p_value)

for domain in ['Cytochrome P450', 'ABC transporter', 'GRAS domain family']:
    D_target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].str.startswith('Dsan'))]['count'].tolist()
    D_medians = np.median(D_target_pfam_values)
    D_q1 = np.percentile(D_target_pfam_values, 10)
    D_q3 = np.percentile(D_target_pfam_values, 90)
    D_stds = np.std(D_target_pfam_values)
    print(D_medians, D_q1, D_q3, D_stds)

    reference_sp = [x for x in phylogeny_class['PAN'] if not x.startswith('Dsan') and not x == 'Zmay' and not x == 'Eoph' and not x.startswith('Pvir')]
    target_pfam_values = plot_pfam_long_df[(plot_pfam_long_df['Pfam']==domain) & (plot_pfam_long_df['Species'].isin(['Ehap','Eory-A','Eory-B','Ecru-A','Ecru-B','Ecru-C']))]['count'].tolist()
    medians = np.median(target_pfam_values)
    q1 = np.percentile(target_pfam_values, 10)
    q3 = np.percentile(target_pfam_values, 90)
    stds = np.std(target_pfam_values)
    print(medians, q1, q3, stds)

    t_stat, p_value = stats.ttest_ind(D_target_pfam_values, target_pfam_values)
    print(t_stat, p_value)