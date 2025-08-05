'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-12-27 15:01:04
LastEditors: Ne0tea
LastEditTime: 2025-06-12 22:59:34
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

# ('C1','C2'),('C1','C3-1'),('C1','C3-2'),('C1','C4'),('C1','C5'),
# ('C2','C3-1'),('C2','C3-2'),('C2','C4'),('C2','C5'),('C3-1','C4'),('C3-1','C5')
pd.set_option('display.max_rows', None)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

Type_color = {'C5-NE':'#264653','C5-S':'#ae2012','C3-1':'#3a5a40','C3-2':'#588157','C4':'#a3b18a',
             'C5-E1':'#94d2bd','C5-SE':'#ee9b00','C5-E2':'#e9d8a6'}
Sp_color = {'C1':'#264653','C2':'#287271','C3-1':'#2a9d8f','C3-2':'#e9c46a','C4':'#f4a261','C5':'#e76f51'}

character_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_all_character.txt'
character_df = pd.read_table(character_file)
#for Cluster paint
# character_df = character_df[(character_df['Cluster'] != 'Outlier') & (character_df['Cluster'] != 'C3-1')]
# character_df = character_df[character_df['Eco_type'].str.startswith('C')]
# print(character_df)
#for Eco_type paint
character_df = character_df[(character_df['Eco_type'] != 'Admix-E1S') &
                            (character_df['Eco_type'] != 'Outlier') &
                            (character_df['Eco_type'] != 'C2') &
                            (character_df['Eco_type'] != 'C1') &
                            (character_df['Cluster'] != 'C3-1') &
                            (character_df['Cluster'] != 'C3-2') &
                            (character_df['Cluster'] != 'C4') &
                            (character_df['Eco_type'] != 'Admix-E2S') &
                            (character_df['Eco_type'] != 'Admix-C4') &
                            (character_df['Eco_type'] != 'Admix-E12')]

character_df['L/W'] = character_df['Leaf length'] / character_df['Leaf width']
character_df['LW'] = character_df['Leaf length'] * character_df['Leaf width']
# character_df = character_df[character_df['Seed weight'] < 100]

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
# print(character_df[classify_column].value_counts())
# print(character_df)
fig = plt.figure(figsize=(8, 10))
# gs = fig.add_gridspec(1, 6,width_ratios=[1.5,1,1,1], height_ratios=[1, 1, 1])
gs = fig.add_gridspec(5, 2)
axs=[]

axs.append(fig.add_subplot(gs[0, 0]))
axs.append(fig.add_subplot(gs[0, 1]))
axs.append(fig.add_subplot(gs[1, 0]))
axs.append(fig.add_subplot(gs[1, 1]))
axs.append(fig.add_subplot(gs[2, 0]))
axs.append(fig.add_subplot(gs[2, 1]))
axs.append(fig.add_subplot(gs[3, 0]))
axs.append(fig.add_subplot(gs[3, 1]))
axs.append(fig.add_subplot(gs[4, 0]))
axs.append(fig.add_subplot(gs[4, 1]))

sns.boxplot(ax=axs[0], y=classify_column, x='Leaf length', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[0], x='Leaf length', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[0], pairs=pairs_dic[classify_column], data=character_df, x='Leaf length', y=classify_column, hue=classify_column, permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[0].grid(True, linestyle='--')
axs[0].set_ylabel('')
axs[0].legend([],[], frameon=False)

sns.boxplot(ax=axs[1], y=classify_column, x='Leaf width', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[1], x='Leaf width', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
print(character_df)
anno = Annotator(axs[1], pairs=pairs_dic[classify_column], data=character_df, x='Leaf width', hue=classify_column, y=classify_column, permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[1].grid(True, linestyle='--')
axs[1].set_ylabel('')
axs[1].legend([],[], frameon=False)

sns.boxplot(ax=axs[2], y=classify_column, x='Stem angle', data=character_df,  hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[2], x='Stem angle', y=classify_column, data=character_df,  hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
# anno = Annotator(axs[2], pairs=pairs_dic[classify_column], data=character_df, x='Stem angle', y=classify_column, permutations=10000,
#                  orient='h', order=list(color_dic[classify_column].keys()))
# anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1)
# anno.apply_and_annotate()
axs[2].grid(True, linestyle='--')
axs[2].set_ylabel('')
axs[2].legend([],[], frameon=False)

sns.boxplot(ax=axs[3], y=classify_column, x='Leaf angle', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[3], x='Leaf angle', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
# anno = Annotator(axs[3], pairs=pairs_dic[classify_column], data=character_df, x='Leaf angle', y=classify_column, permutations=10000,
#                  orient='h', order=list(color_dic[classify_column].keys()))
# anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1)
# anno.apply_and_annotate()
axs[3].grid(True, linestyle='--')
axs[3].set_ylabel('')
axs[3].legend([],[], frameon=False)

sns.boxplot(ax=axs[4], y=classify_column, x='Seed trichome', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[4], x='Seed trichome', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
# anno = Annotator(axs[4], pairs=pairs_dic[classify_column], data=character_df, x='Leaf angle', y=classify_column, permutations=10000,
#                  orient='h', order=list(color_dic[classify_column].keys()))
# anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1)
# anno.apply_and_annotate()
axs[4].grid(True, linestyle='--')
axs[4].set_ylabel('')
axs[4].legend([],[], frameon=False)


sns.boxplot(ax=axs[5], y=classify_column, x='L/W', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[5], x='L/W', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[5], pairs=pairs_dic[classify_column], data=character_df, x='L/W', y=classify_column, hue=classify_column, permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[5].grid(True, linestyle='--')
axs[5].set_ylabel('')

sns.boxplot(ax=axs[6], y=classify_column, x='Seed length', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[6], x='Seed length', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[6], pairs=pairs_dic[classify_column], data=character_df, x='Seed length', hue=classify_column,y=classify_column, permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[6].grid(True, linestyle='--')
axs[6].set_ylabel('')

sns.boxplot(ax=axs[7], y=classify_column, x='Seed width', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[7], x='Seed width', y=classify_column, data=character_df, hue=classify_column,order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[7], pairs=pairs_dic[classify_column], data=character_df, x='Seed width', hue=classify_column,y=classify_column, permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[7].grid(True, linestyle='--')
axs[7].set_ylabel('')

sns.boxplot(ax=axs[8], y=classify_column, x='Seed weight', data=character_df, hue=classify_column,order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[8], x='Seed weight', y=classify_column, data=character_df, hue=classify_column,order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[8], pairs=pairs_dic[classify_column], data=character_df, x='Seed weight', y=classify_column, hue=classify_column,permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[8].grid(True, linestyle='--')
axs[8].set_ylabel('')

sns.boxplot(ax=axs[9], y=classify_column, x='LW', data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
            fill=False,
            palette=color_dic[classify_column], orient='h', showfliers=False)
sns.stripplot(ax=axs[9], x='LW', y=classify_column, data=character_df, hue=classify_column, order=list(color_dic[classify_column].keys()),
              palette=color_dic[classify_column], size=3, jitter=True)
anno = Annotator(axs[9], pairs=pairs_dic[classify_column], data=character_df, x='LW', y=classify_column, hue=classify_column,permutations=700,
                 orient='h', order=list(color_dic[classify_column].keys()), verbose=False)
anno.configure(test='t-test_ind', text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
anno.apply_and_annotate()
axs[9].grid(True, linestyle='--')
axs[9].set_ylabel('')


line_colors = []
for ax in axs:
    for i, patch in enumerate(ax.patches):
        cur_color = patch.get_facecolor()
        line_colors.append(cur_color)
        patch.set_edgecolor(cur_color)
        patch.set_linewidth(1.5)
        patch.set_facecolor('none')

plt.tight_layout()
plt.show()
# plt.savefig(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_character_'+classify_column+'.pdf')