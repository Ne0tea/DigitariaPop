'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-21 15:00:35
LastEditors: Ne0tea
LastEditTime: 2025-03-28 18:37:07
'''
import pandas as pd

def count_values_by_type(df, column_types):
    All_snp_haplotype_df = pd.DataFrame(columns=['Snp','Type','Miss','Ref','Hyp','Hom']) 
    for i in df.columns:
        cur_df=df[i].to_frame()
        cur_df['Type'] = column_types
        for group, gdf in cur_df.groupby('Type'):
            cur_type_dic = gdf[i].value_counts().to_dict()
            cur_dic = {'Snp':i, 'Type':group}
            for var_type in cur_type_dic:
                if var_type == -1:
                    var_str = 'Miss'
                elif var_type == 0:
                    var_str = 'Ref'
                elif var_type == 1:
                    var_str = 'Hyp'
                elif var_type == 2:
                    var_str = 'Hom'
                cur_dic[var_str] = cur_type_dic[var_type]
            cur_df = pd.DataFrame(cur_dic, index=[0])
            All_snp_haplotype_df = pd.concat([All_snp_haplotype_df, cur_df], ignore_index=True)
    return All_snp_haplotype_df

def main(snp_hap, snp_012_pos_file, snp_012_ind_file, herbicide_type_file):
    with open(herbicide_type_file, 'r') as herbf:
        herb_type_dic = {}
        for i in herbf:
            line = i.strip().split('\t')
            herb_type_dic[line[5]] = line[9]
    type_list = []
    valid_index = []
    with open(snp_012_ind_file, 'r') as indf:
        for idx,i in enumerate(indf):
            line = i.strip().split('\t')
            if line[0] not in herb_type_dic: continue
            valid_index.append(idx)
            type_list.append(herb_type_dic[line[0]])
    print(len([ x for x in type_list if x == 'Resistant']), len([ x for x in type_list if x == 'Sensitive']))
    with open(snp_012_pos_file, 'r') as herbf:
        herb_pos_list = []
        for i in herbf:
            line = i.strip().split('\t')
            herb_pos_list.append(str(line[0])+'_'+str(line[1]))
    snp_hap_df = pd.read_table(snp_hap, header=None, sep='\t',index_col=0)
    snp_hap_df = snp_hap_df.loc[valid_index]

    stat_df = count_values_by_type(snp_hap_df, type_list)
    stat_df['Snp'] = stat_df['Snp'].apply(lambda x: herb_pos_list[x-1])
    stat_df.fillna(0, inplace=True)

    classify_stat_df = pd.DataFrame(columns=['Snp','R_Miss','R_Ref','R_Hyp','R_Hom','S_Miss','S_Ref','S_Hyp','S_Hom'])
    for group, i in stat_df.groupby('Snp'):
        R_line = i[i['Type'] == 'Resistant']
        R_line = R_line.rename(columns={'Miss':'R_Miss','Ref':'R_Ref','Hyp':'R_Hyp','Hom':'R_Hom'})
        R_line.drop(columns=['Type'], inplace=True)
        S_line = i[i['Type'] == 'Sensitive']
        S_line = S_line.rename(columns={'Miss':'S_Miss','Ref':'S_Ref','Hyp':'S_Hyp','Hom':'S_Hom'})
        S_line.drop(columns=['Type'], inplace=True)
        cur_line = pd.merge(left=R_line,right=S_line, on='Snp')
        classify_stat_df = pd.concat([classify_stat_df,cur_line], ignore_index=True)
    classify_stat_df.to_csv(snp_hap + r'.stat.txt', sep='\t', index=False)

if __name__ == "__main__":
    # snp_012_hap_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.herbicide.candidategene.012'
    # snp_012_pos_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.herbicide.candidategene.012.pos'
    # snp_012_ind_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.herbicide.candidategene.012.indv'
    # herbicide_type_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_HB_RR.txt'

    snp_012_hap_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_unfilt.012'
    snp_012_pos_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_unfilt.012.pos'
    snp_012_ind_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_unfilt.012.indv'
    herbicide_type_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_HB_RR.txt'

    # snp_012_hap_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_indel_unfilt.012'
    # snp_012_pos_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_indel_unfilt.012.pos'
    # snp_012_ind_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\ALS\DsALS_indel_unfilt.012.indv'
    # herbicide_type_txt = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_HB_RR.txt'

    main(snp_012_hap_file, snp_012_pos_file, snp_012_ind_file, herbicide_type_txt)