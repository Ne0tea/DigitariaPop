'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-02-21 11:30:54
LastEditors: Ne0tea
LastEditTime: 2025-07-18 16:01:36
'''
import pandas as pd
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy import stats
import numpy as np
from matplotlib.patches import Polygon

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'
pd.set_option('display.max_rows', None)
Outerlier_list = ['W-100', '15-18', 'W-9']
Digitaria_platte={'C3-2':'#2a9d8f','C3-1':'#2a9111','C1':'#264653', 'C2':'#287271','C4':'#f4a261','C5':'#e76f51','Outlier':'#000000','Svir':'#000000'}


def draw_sp_radar(data):
    categories = list(data.keys())
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False)

    medians = [np.mean(data[cat]) for cat in categories]
    q1 = [np.percentile(data[cat], 25) for cat in categories]
    q3 = [np.percentile(data[cat], 75) for cat in categories]
    stds = [np.std(data[cat]) for cat in categories]
    print(stds)

    angles = np.concatenate((angles, [angles[0]]))
    medians = np.concatenate((medians, [medians[0]]))
    q1 = np.concatenate((q1, [q1[0]]))
    q3 = np.concatenate((q3, [q3[0]]))
    stds = np.concatenate((stds, [stds[0]]))
    angles = (angles + np.pi/2) % (2 * np.pi)

    fig, ax = plt.subplots(figsize=(4, 4), subplot_kw={'polar': True})
    for level in np.arange(0, 110, 10):
        verts = [(angle, level) for angle in angles]
        poly = Polygon(verts,
                     closed=True,
                     fill=False,
                     edgecolor='#dddddd',
                     linewidth=0.8 if level%20==0 else 0.5,
                     alpha=0.8 if level%20==0 else 0.4,
                     linestyle='-' if level%20==0 else ':')
        ax.add_patch(poly)

    for angle in angles[:-1]:
        ax.plot([angle, angle], [0, 100], color='#aaaaaa', linewidth=1, alpha=0.7)

    ax.plot(angles, medians, color='#2a4d80', linewidth=2.5)

    for angle, median_val, q1_val, q3_val, std_value in zip(angles, medians, q1, q3, stds):
        lower_error = median_val - q1_val
        upper_error = q3_val - median_val

        ax.errorbar(angle, median_val,
                   yerr=[[lower_error], [upper_error]],
                    # yerr=[[std_value], [std_value]],
                   ecolor='#b6d6d6',
                   elinewidth=1.5,
                   capsize=5,
                   capthick=1.5,
                   alpha=0.8)

    ax.set_rlabel_position(90)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=10)
    ax.set_yticks([0,20,40,60])
    ax.set_ylim(0, 60)

    plt.tight_layout()
    # plt.show()
    plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_sp_GR50.pdf', format='pdf')

def add_regression_stats(ax, x, y, color):
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_clean, y_clean)

    stats_text = (
        f"Slope = {slope:.2f}\n"
        f"Std err = {std_err:.2f}\n"
        f"Corr = {r_value:.2f}"
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

material_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Final_material_loction_cluster.txt'
herbicide_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt'
year_color_platte = {2013:'#fb8500', 2015:'#ffb703',2023:'#219ebc'}

herbicide_df = pd.read_csv(herbicide_txt,sep='\t')
material_df = pd.read_csv(material_file,sep='\t')

material_heribicide_df = pd.merge(herbicide_df, material_df, on='ID',)
GR_50_df = material_heribicide_df[material_heribicide_df['Als_GR50'].notna()]
# region filter
# GR_50_df = GR_50_df.loc[(GR_50_df['Longitude']>=108) & (GR_50_df['Longitude']<=120) & (GR_50_df['Latitude']>=32) & (GR_50_df['Latitude']<=40),]
# Ecotype year filter
year3_C5_df = GR_50_df.loc[(GR_50_df['Cluster_in_material']=='C5') & (GR_50_df['Years'].isin([2013,2015,2023])),]
year3_C5_df = year3_C5_df.loc[(year3_C5_df['Longitude']>=108) & (year3_C5_df['Longitude']<=120) & (year3_C5_df['Latitude']>=32) & (year3_C5_df['Latitude']<=40),]
year3_C5_df.loc[:, 'log_GR90'] = year3_C5_df['Als_GR90'].apply(np.log1p)
year3_C5_df.loc[:, 'log_GR50'] = year3_C5_df['Als_GR50'].apply(np.log1p)
# GR_50_df = GR_50_df.loc[GR_50_df['Als_GR90']<1000,]

# # # draw all material year distribution
# pie_data = Counter(material_heribicide_df['Years'])
# sales,channels = pie_data.keys(), pie_data.values()
# print(sum(channels))
# plt.figure(figsize=(4, 3))
# plt.pie(channels, labels=sales, autopct='%1.1f%%', startangle=60, colors=['#e76f51','#f4a261', '#e9c46a', '#264653','#2a9d8f'])

# plt.show()
# # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_Years_distribution.pdf', format='pdf')
# # plt.close()

# # draw all material Cluster distribution
# pie_data = Counter(material_heribicide_df.loc[GR_50_df['Eco_type'].isin(['C1','C2','C3-1','C3-2','C4','C5-E1','C5-E2','C5-NE','C5-S']),'Eco_type'])
# sales,channels = pie_data.keys(), pie_data.values()
# print(sum(channels))
# plt.figure(figsize=(4, 3))
# plt.pie(channels, labels=sales, autopct='%1.1f%%', startangle=60)

# plt.show()
# # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_Eco_distribution.pdf', format='pdf')
# # plt.close()

# # draw C5 3years distribution
# cal_col = 'log_GR90'
# year_color_platte = {2013:'#264653', 2015:'#e9c46a',2016:'#2a9d8f',2021:'#f4a261',2023:'#e76f51'}
# plt.figure(figsize=(3, 5))
# # sns.histplot(data=GR_50_df, x="Als_GR50",hue='Years', hue_order=[2015,2016,2013,2021,2023], bins=60, palette=year_color_platte, multiple='layer', alpha=0.5)
# box_chart = sns.boxplot(data=year3_C5_df, x='Years', y=cal_col, hue='Years', fill=False, 
#                         palette=year_color_platte, 
#                         showfliers=False ,gap=0.2, )

# sns.stripplot(x='Years', y=cal_col, data=year3_C5_df, 
#               palette= {'2013':'#264653', '2015':'#e9c46a','2023':'#e76f51'},
#             legend=False,
#             jitter=True, alpha=0.5)
# # plt.ylim(0,400)
# plt.tight_layout()
# # plt.show()
# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_GR90_C5_boxplot.pdf', format='pdf')
# # plt.close()

# draw Ecotype 3years distribution
cal_col='GR50'
material_heribicide_df['log_GR90'] = material_heribicide_df.loc[:,'Als_GR90'].apply(lambda x: np.log1p(x))
Eco_ord = ['C1','C2','C3-2','C4','C5']
cluster_df = material_heribicide_df[material_heribicide_df['Cluster_y'].isin(Eco_ord)]
cluster_gr50_dic = cluster_df.groupby('Cluster_y')['Als_GR50'].apply(np.array).to_dict()
cluster_df = cluster_df[cluster_df['Als_GR90'] < 1000]
cluster_gr90_dic = cluster_df.groupby('Cluster_y')['Als_GR90'].apply(np.array).to_dict()
categories = list(cluster_gr90_dic.keys())

xx=[]
for cat in ['C2','C4','C5']:
    xx.extend(cluster_gr50_dic[cat])
print(np.median(xx))

for cat in categories:
    mediansxxx = np.median(cluster_gr50_dic[cat])
    print(cat, mediansxxx)
draw_sp_radar(cluster_gr50_dic)

# draw GR50/GR90 log ratio
# year3_C5_df = year3_C5_df.sort_values(by='log_GR90', ascending=False)
# year3_C5_df = year3_C5_df[~year3_C5_df['ID'].isin(Outerlier_list)]
# sns.set_theme(style="darkgrid")

# g = sns.lmplot(
#     x="log_GR50",
#     y="log_GR90",
#     col="Years",
#     hue='Years',
#     data=year3_C5_df,
#     palette = year_color_platte,
#     height=5,
#     aspect=0.6
# )
# g.set_titles("{col_name}")
# for i, ax in enumerate(g.axes.flat):
#     names = ax.get_title()
#     if names == '2013':
#         subset = year3_C5_df[year3_C5_df["Years"] == 2013]
#         color = year_color_platte[2013]
#     elif names == '2015':
#         subset = year3_C5_df[year3_C5_df["Years"] == 2015]
#         color = year_color_platte[2015]
#     elif names == '2023':
#         subset = year3_C5_df[year3_C5_df["Years"] == 2023]
#         color = year_color_platte[2023]
#     else:
#         continue

#     add_regression_stats(ax, subset["log_GR50"], subset["log_GR90"], color)

# # plt.tight_layout()
# plt.show()
# # plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Herbicide_GR_logvs.pdf', format='pdf')