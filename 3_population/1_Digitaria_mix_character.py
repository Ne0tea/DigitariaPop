'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-12-27 15:01:04
LastEditors: Ne0tea
LastEditTime: 2025-04-15 20:57:19
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

# ('C1','C2'),('C1','C3-1'),('C1','C3-2'),('C1','C4'),('C1','C5'),
# ('C2','C3-1'),('C2','C3-2'),('C2','C4'),('C2','C5'),('C3-1','C4'),('C3-1','C5')
pd.set_option('display.max_columns', None)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

Type_color = {'C5-NE':'#264653','C5-S':'#ae2012','C3-1':'#3a5a40','C3-2':'#588157','C4':'#a3b18a',
             'C5-E1':'#94d2bd','C5-SE':'#ee9b00','C5-E2':'#e9d8a6'}
Sp_color = {'C1':'#264653','C2':'#287271','C3-1':'#2a9d8f','C3-2':'#e9c46a','C4':'#f4a261','C5':'#e76f51'}

character_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_all_character.txt'
character_df = pd.read_table(character_file)

character_df['L/W'] = character_df['Seed length'] / character_df['Seed width']
character_df = character_df[character_df['Seed weight'] < 100]
character_df['Tiller angle'] = character_df['Tiller angle'].fillna(0)
character_df['Tiller angle'] = character_df['Tiller angle'].astype(int)

classify_column = 'Class28T'
# classify_column = 'Cluster'
color_dic={'Cluster':{'C1':'#264653','C2':'#287271','C3-2':'#2a9d8f','C4':'#f4a261','C5':'#e76f51'},
           'Eco_type':{'C5-S':'#ae2012','C5-SE':'#ee9b00','C5-E1':'#94d2bd','C5-E2':'#e9d8a6','C5-NE':'#005f73'},
           'Class28T':{'C5-S':'#ae2012','C5-E1':'#94d2bd','C5-E2':'#e9d8a6','C5-NE':'#005f73'}}
pairs_dic={'Cluster':[('C1','C2'),('C2','C3-2'), 
                      ('C3-2','C4'),('C3-2','C5'),
                      ('C4','C5')],
           'Eco_type':[('C5-E1','C5-E2'),('C5-SE','C5-S'),('C5-NE','C5-E2'),('C5-NE','C5-S'),('C5-E1','C5-SE')],
           'Class28T':[('C5-E1','C5-E2'),('C5-NE','C5-E2'),('C5-NE','C5-S'),('C5-E1','C5-S')]}
Eco_admix_order = ['C4','Admix-C4','C5-E1','Admix-E1S','Admix-E12','C5-E2','Admix-E2S','C5-S','C5-NE']
# print(character_df[classify_column].value_counts())
print(character_df)
# fig = plt.figure(figsize=(3, 6))

# sns.boxplot(y=classify_column, x='Stem angle', data=character_df, hue=classify_column, 
#             order=Eco_admix_order,
#             fill=False,
#             # palette=color_dic[classify_column],
#             orient='h', showfliers=False)
# sns.stripplot(x='Stem angle', y=classify_column, data=character_df, hue=classify_column, 
#               order=Eco_admix_order,
#             #   palette=color_dic[classify_column],
#               size=3, jitter=True)

# plt.tight_layout()
# plt.show()

for character in ['Stem color', 'Stem trichome', 'Leaf trichome']:
    grouped = character_df.groupby(['Class28T', character]).size().unstack(fill_value=0)
    percentages = grouped.div(grouped.sum(axis=1), axis=0) * 100

    ax = percentages.plot(kind='bar', stacked=True, figsize=(10, 6))
    plt.title(f'Percentage Composition of {character} by Class28T')
    plt.ylabel('Percentage (%)')
    plt.xlabel('Class28T')

    # 添加百分比标签
    for container in ax.containers:
        ax.bar_label(container, label_type='center', fmt='%.1f%%')

    plt.legend(title=character, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_character_'+classify_column+'.pdf')