'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-07 17:48:33
LastEditors: Ne0tea
LastEditTime: 2025-04-07 22:00:40
'''
import sys
import pandas as pd
signi_file = sys.argv[1]
anno_file = sys.argv[2]
# signi_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.emmax.3.covkin-TopTenPC.assoc.txt.sign7.txt'
signi_df = pd.read_table(signi_file, usecols=[0,1,2], sep='\s+', names=['chr', 'name', 'pos'], engine='python')

signi_df['Chr_name'] = signi_df['chr'].apply(lambda x: 'Chr0' + str(x) if x < 10 else 'Chr' + str(x))
signi_df['pos_region'] = signi_df['pos'].apply(lambda x: range(x-500000, x+500000))
pos_region_list = signi_df['pos_region'].values
with open(anno_file, 'r') as f:
    for i in f:
        line=i.strip().split()
        chr_name, snp_pos = line[0], int(line[1])
        if chr_name in signi_df['Chr_name'].values:
            for pos_region in pos_region_list:
                if snp_pos in pos_region:
                    print(i)
