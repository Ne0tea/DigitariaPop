'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-19 14:08:34
LastEditors: Ne0tea
LastEditTime: 2025-07-26 21:47:08
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
hap_color={'Hap1':'#f5b940','Hap2':'#005d5c'}
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams[ 'font.family'] = 'Arial'

ehh_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\Chr04_54407074.ehh'
# S_ehh_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Chr20_S_34032644.ehh'

R_ehh_df = pd.read_table(ehh_file, sep=' ', skipfooter=2, names=['Snp_name', 'genetic_loc', 'phy_loc','Ref','Alt'])
R_ehh_long_df = R_ehh_df.melt(
    id_vars=[col for col in R_ehh_df.columns if col not in ['Ref', 'Alt']],
    value_vars=['Ref', 'Alt'],
    var_name='Type',
    value_name='Allele'
)
R_ehh_long_df['genetic_loc'] = R_ehh_long_df['genetic_loc'].astype(float)
R_ehh_long_df = R_ehh_long_df[(R_ehh_long_df['genetic_loc'] > 54.3) & (R_ehh_long_df['genetic_loc'] < 54.5)]

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
plt.xlim([54.35, 54.45])
plt.xticks(ticks=[54.35, 54.4, 54.45])

# plt.show()
plt.savefig(ehh_file+".pdf", format='pdf')