'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-21 11:30:54
LastEditors: Ne0tea
LastEditTime: 2025-07-15 20:47:38
'''
import pandas as pd
import seaborn as sns
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
Outerlier_list = ['W-100', '15-18', 'W-9', 'MMT-13-10', 'MMT-13-20']
# year_color_platte_GR90 = {'2013':'#d8f3dc', '2015':'#95d5b2','2016':'#52b788', '2021':'#2d6a4f', '2023':'#1b4332'}
# year_color_platte_GR50 = {'2013':'#ffee99', '2015':'#ffe566','2016':'#ffe14c', '2021':'#ffd819', '2023':'#fdc921'}
year_color_platte_GR90 = {2013:'#d8f3dc', 2015:'#95d5b2', 2016:'#52b788', 2021:'#2d6a4f', 2023:'#1b4332'}
year_color_platte_GR50 = {2013:'#ffe169', 2015:'#fad643', 2016:'#edc531', 2021:'#dbb42c', 2023:'#c9a227'}
year_ord = [2013, 2015, 2016, 2021, 2023]

def draw_boxplot(plot_df):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    # plot_df['group'] = plot_df['Year'].astype(str)
    # plot_df = plot_df[plot_df['Year'].isin(year_ord)]
    # for herbicidetype, color_paltte, ylimset in zip(['GR50', 'GR90'],[year_color_platte_GR50, year_color_platte_GR90],[100, 400]):
    #     plt.figure(figsize=(3, 5))
    #     ax = sns.stripplot(x='Year', y=herbicidetype, hue='Year', data=plot_df,
    #                     #  edgecolor="black",   # 边线颜色
    #                     linewidth=0,
    #                   palette=color_paltte,
    #                 # order=year_ord,
    #                 # legend=False,
    #                 alpha=0.5, edgecolor='black')
    #     plt.ylim(0,ylimset)
    #     plt.show()
    #     # plt.savefig(fr'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\tmp2_{herbicidetype}.pdf', format='pdf')
    #     plt.figure(figsize=(3, 5))
    #     sns.regplot(data=plot_df,x='Year', y=herbicidetype, scatter=False, order=1,
    #             line_kws={'color': '#b6d6d6', 'linestyle': '-'})
    #     plt.ylim(0,ylimset)
    #     plt.grid(True, linestyle='--') 
    #     # plt.xticks(year_ord, year_ord)
    #     # plt.tight_layout()
    #     plt.legend(loc='upper left')
    #     plt.show()
    #     # plt.savefig(fr'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\tmp1_{herbicidetype}.pdf', format='pdf')

    for herbicidetype, color_paltte, ylimset in zip(['GR50', 'GR90'],[year_color_platte_GR50, year_color_platte_GR90],[100, 400]):
        plt.figure(figsize=(3, 3))
        ax = sns.scatterplot(x='Year', y=herbicidetype, hue='Year', data=plot_df,
                        #  edgecolor="black",   # 边线颜色
                        linewidth=0,
                      palette=color_paltte,
                    # order=year_ord,
                    # legend=False,
                    alpha=0.5, edgecolor='black')
        sns.regplot(data=plot_df,x='Year', y=herbicidetype, scatter=False, order=1,
                line_kws={'color': '#b6d6d6', 'linestyle': '-', "linewidth": 2})
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            plot_df["Year"], plot_df[herbicidetype], 'greater'
        )
        r_squared = r_value**2
        text = f"R² = {r_squared:.2f}, p = {p_value:.2e}, " \
               f"F = {(r_squared/(1-r_squared)) * (len(plot_df)-2):.1f}"
        ax.text(
            0.05, 0.95, text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top'
        )
        ax.yaxis.set_major_locator(LinearLocator(3))
        plt.xticks([2013,2018,2023])
        plt.ylim(0,ylimset)
        # plt.grid(True, linestyle='--') 
        plt.legend().set_visible(False)
        plt.show()
        # plt.savefig(fr'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Perinal_{herbicidetype}.pdf', format='pdf')

def main(material_file, herbicide_txt):
    herbicide_df = pd.read_csv(herbicide_txt,sep='\t').set_index('ID')
    material_df = pd.read_csv(material_file,sep='\t')
    herbicide_dic = herbicide_df['Als_GR90'].to_dict()
    herbicide_GR50_dic = herbicide_df['Als_GR50'].to_dict()

    material_df['GR90'] = material_df['ID'].map(herbicide_dic)
    material_df['GR50'] = material_df['ID'].map(herbicide_GR50_dic)
    material_df = material_df[(material_df['Longitude'] > 110) & (material_df['Longitude'] < 125) & \
                              (material_df['Latitude'] > 32) & (material_df['Latitude'] < 40)]
    material_df.dropna(inplace=True)
    material_df = material_df[material_df['Cluster'] == 'C5']
    material_df.reset_index(drop=True,inplace=True)
    material_df = material_df[~material_df['ID'].isin(Outerlier_list)]
    material_df = material_df[abs(material_df['GR90'] - material_df['GR50']) > 5]
    print(material_df.sort_values(by='GR50'))
    # material_df = material_df[((material_df['GR50']>=62) | (material_df['GR50']<58)) & (~material_df['Eco_type'].str.contains('Admix'))].reset_index(drop=True)
    # print(material_df.sort_values(by='GR50'))
    draw_boxplot(material_df)

if __name__ == '__main__':
    ####DROP!
    material_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Final_material_loction_cluster.txt'
    herbicide_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
    main(material_file, herbicide_txt)

