'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-23 23:25:18
LastEditors: Ne0tea
LastEditTime: 2025-07-28 18:00:58
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
year_color_platte = {2013: "#f2c22f", 2015: "#b6d6d6", 2023: "#2d674d"}

def draw_facegrid(plot_df):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    fig = plt.figure(figsize=(5, 2))
    gs = fig.add_gridspec(nrows=2, ncols=3, height_ratios=[8, 1], hspace=0.05)
    axes_hist = []
    for i, chr in enumerate(plot_df["Year"].unique()):
        ax = fig.add_subplot(gs[0, i])
        sns.histplot(
            data=plot_df[plot_df["Year"] == chr],
            x="Freq",
            stat='percent',
            kde=True,
            bins=100,
            color=year_color_platte[chr],
            ax=ax
        )
        ax.set_xlim(0, 0.25)
        ax.set_ylim(0, 50)
        ax.set_title(f"{chr}")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.tick_params(axis='y', which='major', pad=1)
        axes_hist.append(ax)

    axes_strip = []
    for i, chr in enumerate(plot_df["Year"].unique()):
        ax = fig.add_subplot(gs[1, i])
        sns.boxplot(
            data=plot_df[plot_df["Year"] == chr],
            x="Freq",
            # y="Freq",
            orient='h',
            color=year_color_platte[chr],
            fill=False, 
            showfliers=False,
            ax=ax
        )
        ax.set_xlabel("Freq")
        axes_strip.append(ax)

    for ax_hist, ax_strip in zip(axes_hist, axes_strip):
        xmin, xmax = ax_hist.get_xlim()
        ax_strip.set_xlim(xmin, xmax)

    plt.suptitle("", y=0.95)
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_allele_freq.pdf', format='pdf')
    # plt.show()

def draw_propotion(plot_df):
    # 设置图表样式
    sns.set_theme(style="darkgrid")
    g = sns.displot(
        x="Freq",
        row="Year",
        hue='Year',
        data=plot_df,
        stat='percent',
        palette = year_color_platte,
        bins=80,
        height=7,
        facet_kws={"sharey": False},
        aspect=0.8
    )
    g.set_titles("{col_name}")
    plt.xlim(0, 0.25)
    plt.xlabel("Freq", fontsize=12)
    plt.ylabel("Allele Frequency", fontsize=12)

    # plt.savefig("af_boxplot.png", dpi=300, bbox_inches="tight")
    plt.show()

def main(Freq_2013_file, Freq_2015_file, Freq_2023_file):
    Freq_2013_df = pd.read_csv(Freq_2013_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_2013_df['Year'] = 2013
    Freq_2015_df = pd.read_csv(Freq_2015_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_2015_df['Year'] = 2015
    Freq_2023_df = pd.read_csv(Freq_2023_file, sep="\t", names=['Chr', 'Loc', 'Ref', 'Alt', 'Freq','Info'])
    Freq_2023_df['Year'] = 2023

    df_combined = pd.concat([Freq_2013_df, Freq_2015_df, Freq_2023_df])
    # df_combined = pd.concat([Freq_2013_df, Freq_2023_df])
    df_cleaned_combined = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                   .filter(lambda x: x['Freq'].ne(0).any())

    if 0:
        # if for stat allele freq
        Cleaned_new_in13 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2013]['Freq'].ne(0)) & \
                               (x[x['Year']==2023]['Freq'].eq(0))
                               )
        Cleaned_new_in23 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2023]['Freq'].ne(0)) & \
                               (x[x['Year']==2013]['Freq'].eq(0))
                               )
        Higher_in23 = df_combined.groupby(['Chr', 'Loc'], group_keys=False) \
                       .filter(lambda x: (x[x['Year']==2023]['Freq'] > x[x['Year']==2013]['Freq']) & \
                               (x[x['Year']==2013]['Freq'].ne(0))
                               )
        Cleaned_new_in23.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_E1_new_in23.csv', sep='\t', index=False)
        Higher_in23.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_2015_2023_E1_Higher_in23.csv', sep='\t', index=False)
        print(Cleaned_new_in23)
        print(Cleaned_new_in13)
        print(Higher_in23)

    df_cleaned_combined['Type'] = df_cleaned_combined['Info'].apply(lambda x: x.split('|')[1])
    # draw_propotion(df_cleaned_combined)
    draw_facegrid(df_cleaned_combined)
if __name__ == "__main__":
    # Freq_2013_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'
    # Freq_2015_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2015_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'
    # Freq_2023_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2023_Ecotype_E1.txt.candidate_SNPs.MsModi.reAF.pos.freq'

    Freq_2013_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2013_Ecotype.candidate_SNPs.MsModi.reAF.pos.freq'
    Freq_2015_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2015_Ecotype.candidate_SNPs.MsModi.reAF.pos.freq'
    Freq_2023_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Perinal_analysis\Perinal_2023_Ecotype.candidate_SNPs.MsModi.reAF.pos.freq'

    main(Freq_2013_file, Freq_2015_file, Freq_2023_file)