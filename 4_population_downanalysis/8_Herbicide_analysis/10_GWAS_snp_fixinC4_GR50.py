'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-23 23:25:18
LastEditors: Ne0tea
LastEditTime: 2025-05-14 15:03:17
'''
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
import seaborn as sns

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
colors = ["#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226"]
columns = [0, 1, 2, 3, 4, 5]
year_color_platte = {2013:'#fb8500', 2015:'#ffb703',2023:'#219ebc'}
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
def draw_scatterplot(plot_df):
    fig = plt.figure(figsize=(9, 3))
    gs = fig.add_gridspec(nrows=1, ncols=3)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    sns.scatterplot(ax=ax0,data=plot_df,x='C4_fix_snp',y='GR50',hue='Ecotype',s=25,
        # color='#023e8a',
        palette=Eco_type_color)
    sns.scatterplot(ax=ax1,data=plot_df,x='C5_fix_snp',y='GR50',hue='Ecotype',s=25,
        # color='#023e8a',
        palette=Eco_type_color)
    sns.scatterplot(ax=ax2,data=plot_df,x='C5_denovo_snp',y='GR50',hue='Ecotype',s=25,
        # color='#023e8a',
        palette=Eco_type_color)
    # plt.ylim(0, 600)
    plt.title('')
    plt.show()


def calculate_loc_fix(subset_df, material_class_df, meta):
    C4_material_list = material_class_df[material_class_df['Ecotype']=='C4'].index.tolist()
    C4_df = subset_df.loc[C4_material_list]
    C4_fix_genotypes = C4_df.mode()
    C4_fix_genotypes_counts = C4_df.value_counts().max()
    if (len(C4_fix_genotypes) > 1) | (C4_fix_genotypes_counts < (len(C4_df) / 2)) | (C4_fix_genotypes[0] == './.'):
        return False

    C5_list = material_class_df[(material_class_df['Class']=='C5')].index.tolist()
    C5_df: pd.Series = subset_df.loc[C5_list]
    C5_fix_genotypes = C5_df.mode()
    C5_fix_genotypes_counts = C5_df.value_counts().max()

    C5_GR50_list = material_class_df[(material_class_df['Class']=='C5') & (~material_class_df['GR50'].isna())].index.tolist()
    # C5_GR50_list = material_class_df[(material_class_df['Class']=='C5')].index.tolist()
    C5_GR50_df: pd.Series = subset_df.loc[C5_GR50_list]

    for idx, value in C5_GR50_df.items():
        if value == C4_fix_genotypes[0] and value != C5_fix_genotypes[0] and value != './.':
            # if idx == '23-2-10_1':
            #     print(idx, meta[['CHROM','POS']])
            C5_GR50_C4_fix_dic[idx] += 1
        elif value == C5_fix_genotypes[0] and value != C4_fix_genotypes[0] and value != './.':
            C5_GR50_C5_fix_dic[idx] += 1
        elif value != C5_fix_genotypes[0] and value != C4_fix_genotypes[0] and value != './.':
            C5_GR50_denovo_dic[idx] += 1

def main(vcf_file, classify_file):
    with open(vcf_file, 'r') as f:
        lines = [line for line in f if not line.lstrip().startswith('##')]
    vcf_df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=0)
    vcf_df= vcf_df.rename(columns={'#CHROM':'CHROM'})
    cols_to_modify = vcf_df.columns[9:]
    vcf_df[cols_to_modify] = vcf_df[cols_to_modify].apply(
        lambda col: col.str.split(':').str[0]
    )

    material_class_df = pd.read_table(classify_file, header=0, sep='\t', index_col=0)
    C5_GR50_list = material_class_df[(material_class_df['Class']=='C5') & (~material_class_df['GR50'].isna())].index.tolist()
    # C5_GR50_list = material_class_df[(material_class_df['Class']=='C5')].index.tolist()

    global C5_GR50_C4_fix_dic
    C5_GR50_C4_fix_dic = {x:0 for x in C5_GR50_list}
    global C5_GR50_C5_fix_dic
    C5_GR50_C5_fix_dic = {x:0 for x in C5_GR50_list}
    global C5_GR50_denovo_dic
    C5_GR50_denovo_dic = {x:0 for x in C5_GR50_list}


    for idx, row_df in vcf_df.iterrows():
        material_df = row_df[cols_to_modify]
        meta = row_df[vcf_df.columns[:9]]
        calculate_loc_fix(material_df, material_class_df, meta)

    C5_GR50_C4_fix_df = pd.DataFrame.from_dict(C5_GR50_C4_fix_dic, orient='index')
    C5_GR50_C4_fix_df.columns = ['C4_fix_snp']
    C5_GR50_C5_fix_df = pd.DataFrame.from_dict(C5_GR50_C5_fix_dic, orient='index')
    C5_GR50_C5_fix_df.columns = ['C5_fix_snp']
    C5_GR50_denovo_df = pd.DataFrame.from_dict(C5_GR50_denovo_dic, orient='index')
    C5_GR50_denovo_df.columns = ['C5_denovo_snp']

    plot_df = pd.merge(material_class_df, C5_GR50_C4_fix_df, left_index=True, right_index=True, how='inner')
    plot_df = pd.merge(plot_df, C5_GR50_C5_fix_df, left_index=True, right_index=True, how='inner')
    plot_df = pd.merge(plot_df, C5_GR50_denovo_df, left_index=True, right_index=True, how='inner')
    # plot_df = plot_df[~plot_df['Ecotype'].str.contains('Admix')]
    print(plot_df.sort_values('GR90', ascending=False))
    # draw_scatterplot(plot_df)

    # plot_df = pd.merge(C5_GR50_C4_fix_df, C5_GR50_C5_fix_df, left_index=True, right_index=True, how='inner')
    # plot_df = pd.merge(plot_df, C5_GR50_denovo_df, left_index=True, right_index=True, how='inner')
    # # plot_df = plot_df[~plot_df['Ecotype'].str.contains('Admix')]
    # print(plot_df.sort_values('C4_fix_snp', ascending=False))
if __name__ == "__main__":
    vcf_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_GR50_uniq.gene.mis.eff.vcf'
    classify_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Material_Ecotype.list'

    main(vcf_file, classify_file)