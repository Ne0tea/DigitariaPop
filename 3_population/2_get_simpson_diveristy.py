'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-06-02 16:20:46
LastEditors: Ne0tea
LastEditTime: 2025-06-17 20:32:30
'''

import pandas as pd
from skbio.diversity import alpha_diversity
from skbio.diversity import get_alpha_diversity_metrics

Digitaria_all_loc_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_all_sample_loc.txt'
Digitaria_all_loc_df = pd.read_table(Digitaria_all_loc_file)
Digitaria_all_loc_df = Digitaria_all_loc_df[Digitaria_all_loc_df['Cluster'] != 'Outlier']
Digitaria_all_loc_df = Digitaria_all_loc_df[Digitaria_all_loc_df['Cluster'] != 'C3-1']
Digitaria_grouped_count = Digitaria_all_loc_df.groupby("Provience")["Cluster"].value_counts().unstack(fill_value=0)

provience_list = Digitaria_grouped_count.index.tolist()  # Provience 列表
cluster_counts = Digitaria_grouped_count.values.tolist()

print(get_alpha_diversity_metrics())
diversity = alpha_diversity('shannon', cluster_counts, ids=provience_list)
# diversity = alpha_diversity('simpson', cluster_counts, ids=provience_list)

Digitaria_grouped_count = pd.concat([Digitaria_grouped_count, diversity], axis=1)
Digitaria_grouped_count['C5_propotion'] = Digitaria_grouped_count['C5'] / Digitaria_grouped_count.sum(axis=1)
print(Digitaria_grouped_count)
print(len(Digitaria_grouped_count[Digitaria_grouped_count.sum(axis=1) > 5]))
print(len(Digitaria_grouped_count[Digitaria_grouped_count['C5_propotion'] > 0.8]))
Digitaria_grouped_count.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Character_compare\Digitaria_provience_Shannon_diversity.csv', index=True)
