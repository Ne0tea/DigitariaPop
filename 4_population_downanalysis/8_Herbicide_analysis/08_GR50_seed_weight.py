'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-21 11:30:54
LastEditors: Ne0tea
LastEditTime: 2025-04-10 16:07:48
'''
import pandas as pd
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
Outerlier_list = ['W-100', '15-18', 'W-9']
year_color_platte = {2013:'#fb8500', 2015:'#ffb703',2023:'#219ebc'}
Eco_type_color={'C5-NE':'#005f73','C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6','C4':'#ffb703',
             'Admix-C4':'#3d405b','Admix-E12':'#81b29a','Admix-E2S':'#f4f1de','Admix-E1S':'#e07a5f',
             'OUT':'#000000'}
def add_regression_stats(ax, x, y, color):
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_clean, y_clean)

    stats_text = (
        f"Slope = {slope:.5f}\n"
        f"Std err = {std_err:.5f}\n"
        f"Corr = {r_value:.5f}"
    )

    ax.text(
        0.95, 0.15,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        color=color,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none')  # 背景框
    )

def draw_scatter_plot(year3_C5_df, columnx, columny):
    sns.set_theme(style="darkgrid")
    # plt.figure(figsize=(5, 5))
    # sns.scatterplot(x='log_GR50', y='log_GR90', hue='Years' ,data=year3_C5_df,
    #             #   palette=year_color_platte,
    #             alpha=0.5)
    g = sns.lmplot(
        x=columnx,
        y=columny,
        col="Class28T",
        hue='Class28T',
        data=year3_C5_df,
        palette = Eco_type_color,
        height=5,
        aspect=0.6
    )
    g.set_titles("{col_name}")
    for i, ax in enumerate(g.axes.flat):
        names = ax.get_title()
        if names == 'C5-NE':
            subset = year3_C5_df[year3_C5_df["Class28T"] == 'C5-NE']
            color = Eco_type_color['C5-NE']
        elif names == 'C5-S':
            subset = year3_C5_df[year3_C5_df["Class28T"] == 'C5-S']
            color = Eco_type_color['C5-S']
        elif names == 'C5-E1':
            subset = year3_C5_df[year3_C5_df["Class28T"] == 'C5-E1']
            color = Eco_type_color['C5-E1']
        elif names == 'C5-E2':
            subset = year3_C5_df[year3_C5_df["Class28T"] == 'C5-E2']
            color = Eco_type_color['C5-E2']
        elif names == 'C4':
            subset = year3_C5_df[year3_C5_df["Class28T"] == 'C4']
            color = Eco_type_color['C4']
        else:
            continue

        add_regression_stats(ax, subset[columnx], subset[columny], color)

    plt.tight_layout()
    plt.show()
    # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_GR_logvs.pdf', format='pdf')


def main(material_file, material_character_file, herbicide_txt):
    herbicide_df = pd.read_csv(herbicide_txt,sep='\t')
    material_df = pd.read_csv(material_file,sep='\t')
    character_df = pd.read_csv(material_character_file,sep='\t')
    material_heribicide_df = pd.merge(herbicide_df, material_df, on='ID',)
    material_GR50_character_df = pd.merge(material_heribicide_df, character_df, on='ID',)[['ID', 'Cluster_x', 'Class28T', 'Years', 'Als_GR50', 'Seed weight', 'Leaf length', 'Leaf width']]
    # material_GR50_character_df = material_GR50_character_df[(material_GR50_character_df['Cluster_x'] == 'C5') & (material_GR50_character_df['Years'].isin([2013,2015,2023]))]
    material_GR50_character_df = material_GR50_character_df[(material_GR50_character_df['Cluster_x'].isin(['C4','C5'])) & (material_GR50_character_df['Class28T'].isin(['C5-NE','C5-S','C5-E1','C5-E2','C4']))]
    # print(material_GR50_character_df)
    # material_GR50_character_df['Als_GR50'] = material_GR50_character_df['Als_GR50'].apply(lambda x: np.log(x))
    for i in ['Seed weight', 'Leaf length', 'Leaf width']:
        draw_scatter_plot(material_GR50_character_df, 'Als_GR50', i)

if __name__ == '__main__':
    material_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Final_material_loction_cluster.txt'
    material_character_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_all_character.txt'
    herbicide_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
    main(material_file, material_character_file, herbicide_txt)
