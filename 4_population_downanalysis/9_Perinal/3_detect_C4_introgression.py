'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-23 23:25:18
LastEditors: Ne0tea
LastEditTime: 2025-04-18 21:54:56
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'

colors = ["#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226"]
columns = [0, 1, 2, 3, 4, 5]
year_color_platte = {2013:'#fb8500', 2015:'#ffb703',2023:'#219ebc'}

def main(Freq_2013_file, Freq_2015_file, Freq_2023_file, Freq_C4_file):
    Freq_2013_df = pd.read_csv(Freq_2013_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_2013_df['Year'] = 2013
    # Freq_2015_df = pd.read_csv(Freq_2015_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    # Freq_2015_df['Year'] = 2015
    Freq_2023_df = pd.read_csv(Freq_2023_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_2023_df['Year'] = 2023
    Freq_C4_df = pd.read_csv(Freq_C4_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_C4_df['Year'] = 'C4'


    # df_combined = pd.concat([Freq_2013_df, Freq_2015_df, Freq_2023_df])
    df_combined = pd.concat([Freq_2013_df, Freq_2023_df, Freq_C4_df])
    df_cleaned_combined = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                   .filter(lambda x: x['Freq'].ne(0).any())

    if 1:
        # if for stat allele freq
        Cleaned_new_in13 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2013]['Freq'].ne(0)) & \
                               (x[x['Year']==2023]['Freq'].eq(0))
                               )
        Cleaned_new_in13.sort_values(['Chr','Loc'], inplace=True)
        Cleaned_new_in23 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2023]['Freq'].ne(0)) & \
                               (x[x['Year']==2013]['Freq'].eq(0))
                               )
        Cleaned_new_in23.sort_values(['Chr','Loc'], inplace=True)
        Higher_in23 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2023]['Freq'] > x[x['Year']==2013]['Freq']) & \
                               (x[x['Year']==2013]['Freq'].ne(0))
                               )
        Higher_in23.sort_values(['Chr','Loc'], inplace=True)
        Intro_from_C4_in23= df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2023]['Freq'] > x[x['Year']==2013]['Freq']) & \
                               (x[x['Year']==2013]['Freq'].eq(0)) & \
                               (x[x['Year']=='C4']['Freq'].ne(0))
                               )
        Intro_from_C4_in23.sort_values(['Chr','Loc'], inplace=True)
        # Cleaned_new_in23.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_E1_new_in23.csv', sep='\t', index=False)
        # Higher_in23.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_E1_Higher_in23.csv', sep='\t', index=False)
        print(Cleaned_new_in23)
        print(Cleaned_new_in13)
        print(Higher_in23)
        print(Intro_from_C4_in23)
        Intro_from_C4_in23.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_E1_Intro_from_C4_in23.csv', sep=',', index=False)
    df_cleaned_combined['Type'] = df_cleaned_combined['Info'].apply(lambda x: x.split('|')[1])

if __name__ == "__main__":
    Freq_2013_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'
    Freq_2015_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2015_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'
    Freq_2023_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2023_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'
    Freq_C4_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Sample_C4.candidate_SNPs.MsModi.reAF.pos.freq'

    main(Freq_2013_file, Freq_2015_file, Freq_2023_file, Freq_C4_file)