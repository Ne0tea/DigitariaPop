'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-07 21:03:32
LastEditors: Ne0tea
LastEditTime: 2025-03-08 16:22:02
'''
from sklearn.cluster import KMeans
from scipy.stats import shapiro, boxcox
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'

Herbicide_resistance_df = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_herbicide_character.txt')
sensitive_sample = Herbicide_resistance_df[Herbicide_resistance_df['Als_GR90']<50]
resistance_sample = Herbicide_resistance_df[(Herbicide_resistance_df['Als_GR90']>=50)]

mean_sensitive = sensitive_sample[sensitive_sample['Cluster_in_material']=='C5']['Als_GR90'].mean()
print(mean_sensitive)
resistance_sample['RI'] = resistance_sample['Als_GR90'] / mean_sensitive
stats = resistance_sample['RI'].describe()

data = np.log1p(resistance_sample['RI']).values.reshape(-1, )

# plt.figure(figsize=(2, 5))
# # sns.histplot(data=GR_50_df, x="Als_GR50",hue='Years', hue_order=[2015,2016,2013,2021,2023], bins=60, palette=year_color_platte, multiple='layer', alpha=0.5)
# box_chart = sns.boxplot(data=resistance_sample, x='Years', y="RI", hue='Years', fill=False, 
#                         # palette=year_color_platte,
#                         showfliers=False ,gap=0.2)

# sns.stripplot(x='Years', y='RI', hue='Years', data=resistance_sample, 
#             #   palette=year_color_platte,
#             legend=False,
#             dodge=True, jitter=True, alpha=0.5, edgecolor='black')
# plt.ylim(0, 60)
# plt.tight_layout()
# plt.show()

# mu = np.mean(data)
# sigma = np.std(data)

# z_scores = (data - mu) / sigma

# low_group = data[z_scores < -1.28]
# high_group = data[z_scores > 1.28]
# mid_group = data[(z_scores >= -1.28) & (z_scores <= 1.28)]
# print(len(low_group), len(high_group), len(mid_group))

# plot_df = pd.DataFrame({
#     'Als_GR50': np.log(Herbicide_resistance_df['Als_GR50']),
#     'Als_GR90': np.log(Herbicide_resistance_df['Als_GR90'])
# })

# g = sns.JointGrid(data=plot_df, x="Als_GR50", y="Als_GR90", height=8)
# g.plot_joint(sns.scatterplot, alpha=0.6, color="#4ECDC4")
# g.plot_marginals(sns.histplot, kde=False, color="#4ECDC4", edgecolor="white", bins=20)
# g.ax_joint.set(xlabel="Als_GR50", ylabel="Als_GR90", title="")
# plt.tight_layout()
# plt.show()