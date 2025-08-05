'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-19 14:08:34
LastEditors: Ne0tea
LastEditTime: 2025-07-17 16:36:16
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
hap_color={'Hap1':'#f5b940','Hap2':'#005d5c'}
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'

# R_ehh_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Chr20_R_34032644.ehh'
S_ehh_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Chr20_S_34032644.ehh'

R_ehh_df = pd.read_table(S_ehh_file, sep=' ', skipfooter=2, names=['Snp_name', 'genetic_loc', 'phy_loc','Ref','Alt'])
R_ehh_long_df = R_ehh_df.melt(
    id_vars=[col for col in R_ehh_df.columns if col not in ['Ref', 'Alt']],  # 保留其他所有列
    value_vars=['Ref', 'Alt'],  # 要转换的列
    var_name='Type',  # 新列名
    value_name='Allele'  # 值列名
)
R_ehh_long_df['genetic_loc'] = R_ehh_long_df['genetic_loc'].astype(float)
R_ehh_long_df = R_ehh_long_df[(R_ehh_long_df['genetic_loc'] > 33.9) & (R_ehh_long_df['genetic_loc'] < 34.1)]


custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

sns.lineplot(
    data=R_ehh_long_df,
    x="genetic_loc", 
    y="Allele", 
    hue="Type",
    # style="sex",  # 可选：不同组使用不同线型
    markers=True,  # 可选：显示数据点标记
    palette={"Ref": "#f5b940", "Alt": "#005d5c"}  # 自定义颜色
)
plt.xticks(ticks=[33.9,34.0,34.1])

plt.show()
# plt.savefig(R_ehh_file+".pdf", format='pdf')