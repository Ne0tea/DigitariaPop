'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-21 11:30:54
LastEditors: Ne0tea
LastEditTime: 2025-04-17 10:51:48
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
Outerlier_list = ['W-100', '15-18', 'W-9']
year_color_platte = {2013:'#fb8500', 2015:'#ffb703',2023:'#219ebc'}
year_ord = [2013, 2015, 2023]

def draw_scatterplot(plot_df):
    sns.scatterplot(
        data=plot_df,
        x='GCount',
        y='Als_GR50',
        s=25,
        color='#023e8a'
    )

    plt.title('Sample vs Gene Count (Seaborn)')
    plt.show()

def main(herbicide_file, Rsistant_sample_hap_file):
    herbicide_df = pd.read_table(herbicide_file,sep='\t').set_index('ID')
    herbicide_df = herbicide_df[(herbicide_df['Cluster_in_material'] == 'C5') &  (~herbicide_df['Eco_cluster'].str.contains('C4'))]
    resistant_hap_df = pd.read_table(Rsistant_sample_hap_file,sep='\t',names=['Vcf_id', 'Hap', 'Gene'])
    sample_resistant_count_df = resistant_hap_df.groupby('Vcf_id').agg(Gene_value=('Gene', list),
                                           GCount=('Gene', 'nunique')
                                          ).reset_index()
    Sample_GR50_hap_df = pd.merge(herbicide_df, sample_resistant_count_df, how='left', on='Vcf_id')[['Vcf_id','Gene_value','GCount','Als_GR50']]
    Sample_GR50_hap_df.fillna(0, inplace=True)
    # print(Sample_GR50_hap_df)
    draw_scatterplot(Sample_GR50_hap_df)

    resistant_hap_df = pd.merge(resistant_hap_df, herbicide_df, how='left', on='Vcf_id')
    resistant_hap_df.dropna(inplace=True)
    gene_resistant_count_df = resistant_hap_df.groupby('Gene').agg(GR50=('Als_GR50', list),
                                           SampleCount=('Vcf_id', ','.join)
                                          ).reset_index()
    gene_resistant_count_df['Mean_GR50'] = gene_resistant_count_df['GR50'].apply(lambda x: np.mean(x))
    gene_resistant_count_df['Median_GR50'] = gene_resistant_count_df['GR50'].apply(lambda x: np.median(x))
    print(gene_resistant_count_df.sort_values('SampleCount', ascending=False))
    # gene_resistant_count_df = resistant_hap_df.groupby('Gene').agg(Gene_value=('Gene', list),
    #                                        GCount=('Gene', 'nunique')
    #                                       ).reset_index()


if __name__ == '__main__':
    herbicide_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
    Rsistant_sample_hap_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\GWAS_GR50_fdgene_Rhap.list'
    main(herbicide_file, Rsistant_sample_hap_file)

