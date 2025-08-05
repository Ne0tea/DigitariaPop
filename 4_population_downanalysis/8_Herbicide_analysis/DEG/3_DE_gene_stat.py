'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-16 23:16:57
LastEditors: Ne0tea
LastEditTime: 2025-07-28 12:29:14
'''
import pandas as pd
from scipy import stats
import numpy as np
# pd.set_option('display.max_rows', None)

def read_stringtie_out(file):
    result_df = pd.read_csv(file, header=0)
    result_df = result_df.set_index('gene_id')
    return result_df

def make_treat_dic(treat_file):
    treat_name_SRR={}
    CK_list = []
    with open(treat_file, 'r') as gf:
        for i in gf:
            line=i.strip().split()
            srr_id = line[0]
            libary_type = line[1].split('|')[1]
            if libary_type == 'CK':
                CK_list.append(srr_id)
                continue
            treat_name = line[1].split('|')[0].split('_')[1]
            if srr_id == '15-2-0-3': continue
            if treat_name in treat_name_SRR:
                treat_name_SRR[treat_name].append(srr_id)
            else:
                treat_name_SRR[treat_name]=[srr_id]

    return treat_name_SRR, CK_list

def get_treatment_expr(raw_expr_df, all_treat_name, CK_list, prefix):
    temp_Sen_id = [x for x in all_treat_name if x.startswith(prefix)]

    cur_count_df = raw_expr_df[temp_Sen_id]
    cur_count_df = cur_count_df.reset_index()
    cur_count_df_long = cur_count_df.melt(id_vars="gene_id",var_name="sample",value_name="FPKM")
    cur_count_df_long[["Var", "Ht", "Time", "Rep"]] = cur_count_df_long["sample"].str.split("-", expand=True)[[0, 2,3,4]]

    cur_CK_df = raw_expr_df[[x for x in CK_list if x.startswith(prefix)]]
    cur_CK_df = cur_CK_df.reset_index()
    cur_CK_df_long = cur_CK_df.melt(id_vars="gene_id",var_name="sample",value_name="FPKM")
    cur_CK_df_long[["Var", "Time", "Rep"]] = cur_CK_df_long["sample"].str.split("-", expand=True)[[0, 2,3]]
    cur_CK_df_long['Ht'] = cur_count_df_long["Ht"].to_list()[0]
    cur_CK_df_long['Time'] = '0h'
    cur_count_df_long = pd.concat([cur_count_df_long, cur_CK_df_long], axis=0)

    return cur_count_df_long


def main(treat_file, DE_gene_file, R_count_file, S_count_file):
    R_count_df = read_stringtie_out(R_count_file)
    S_count_df = read_stringtie_out(S_count_file)
    S_count_df = S_count_df.drop(columns=['15-2-0-3'])

    treat_name_SRR, CK_list = make_treat_dic(treat_file)
    # CK_list.remove('15-2-0-3')

    DE_expr_df = pd.read_table(DE_gene_file, header=0)
    DE_expr_df[['Var', 'Ht', 'Time']] = DE_expr_df["Proj"].str.split("_", expand=True)[[0,1,2]]

    for Ht, Ht_df in DE_expr_df.groupby("Ht"):
        if Ht != 'Y': continue
        cur_R_long_df = get_treatment_expr(R_count_df, treat_name_SRR[Ht], CK_list, '21-17')
        cur_R_long_FPKM_dict = cur_R_long_df.groupby(['gene_id', 'Time'])['FPKM'].apply(list).to_dict()
        cur_S_long_df = get_treatment_expr(S_count_df, treat_name_SRR[Ht], CK_list, '15-2')
        cur_S_long_FPKM_dict = cur_S_long_df.groupby(['gene_id', 'Time'])['FPKM'].apply(list).to_dict()
        DE_in_RS = []
        for gene_id in Ht_df['gene_id'].unique():
        # for gene_id in ['Chr01.1001']:
            cur_R_0h_expr, cur_R_6h_expr, cur_R_24h_expr = cur_R_long_FPKM_dict[(gene_id, '0h')], cur_R_long_FPKM_dict[(gene_id, '6h')], cur_R_long_FPKM_dict[(gene_id, '24h')]
            cur_S_0h_expr, cur_S_6h_expr, cur_S_24h_expr = cur_S_long_FPKM_dict[(gene_id, '0h')], cur_S_long_FPKM_dict[(gene_id, '6h')], cur_S_long_FPKM_dict[(gene_id, '24h')]
            _, cur_0h_p = stats.ttest_ind(cur_R_0h_expr, cur_S_0h_expr, alternative='greater')
            _, cur_6h_p = stats.ttest_ind(cur_R_6h_expr, cur_S_6h_expr, alternative='greater')
            # _, FPKM_6h_p = stats.ttest_ind(cur_R_6h_expr, cur_R_0h_expr, alternative='greater')
            # _, FPKM_24h_p = stats.ttest_ind(cur_R_24h_expr, cur_R_0h_expr, alternative='greater')
            _, cur_24h_p = stats.ttest_ind(cur_R_24h_expr, cur_S_24h_expr, alternative='greater')
            print(gene_id, cur_6h_p, cur_24h_p)
            if np.mean(cur_R_6h_expr + cur_S_6h_expr) > 5 and (cur_6h_p < 0.05 or cur_0h_p < 0.05 or cur_24h_p < 0.05):
                DE_in_RS.append(gene_id)

            # if np.mean(cur_R_24h_expr + cur_S_24h_expr) > 10 and cur_24h_p < 0.05:
            #     DE_in_RS.append(gene_id)
        creat_df = DE_expr_df[DE_expr_df['gene_id'].isin(DE_in_RS) & (DE_expr_df['Ht'] == Ht)]
        creat_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_RvS_alltime_dif_gene.csv',index=False)

        # Ht_df['log2BaseMean'] = np.log2(Ht_df['baseMean'])
        # DE_gene_list = []
        # for gene, gene_df in Ht_df.groupby("gene_id"):
        #     Res_24_FC=gene_df.loc[(gene_df['Time']=='24h') & (gene_df['Var']=='Resistance')]['log2FoldChange'].values[0]
        #     Res_6_FC=gene_df.loc[(gene_df['Time']=='6h') & (gene_df['Var']=='Resistance')]['log2FoldChange'].values[0]
        #     Sen_24_FC=gene_df.loc[(gene_df['Time']=='24h') & (gene_df['Var']=='Sensitive')]['log2FoldChange'].values[0]
        #     Sen_6_FC=gene_df.loc[(gene_df['Time']=='6h') & (gene_df['Var']=='Sensitive')]['log2FoldChange'].values[0]
        #     # Res_24_FPKM=cur_R_long_df.loc[(cur_R_long_df['gene_id']==gene) & (cur_R_long_df['Time']=='24h'),'FPKM'].mean()
        #     Res_24_FPKM = np.mean(cur_R_long_FPKM_dict[(gene, '24h')])
        #     # Res_6_FPKM=cur_R_long_df.loc[(cur_R_long_df['gene_id']==gene) & (cur_R_long_df['Time']=='6h'),'FPKM'].mean()
        #     Res_6_FPKM = np.mean(cur_R_long_FPKM_dict[(gene, '6h')])
        #     # Sen_24_FPKM=cur_S_long_df.loc[(cur_S_long_df['gene_id']==gene) & (cur_S_long_df['Time']=='24h'),'FPKM'].mean()
        #     Sen_24_FPKM = np.mean(cur_S_long_FPKM_dict[(gene, '24h')])
        #     # Sen_6_FPKM=cur_S_long_df.loc[(cur_S_long_df['gene_id']==gene) & (cur_S_long_df['Time']=='6h'),'FPKM'].mean()
        #     Sen_6_FPKM = np.mean(cur_S_long_FPKM_dict[(gene, '6h')])
        #     if max([Res_24_FPKM, Res_6_FPKM, Sen_24_FPKM, Sen_6_FPKM]) < 10:
        #         continue
        #     print(gene)
        #     if Res_24_FPKM > Sen_24_FPKM and Res_6_FPKM > Sen_6_FPKM and Res_6_FC > 1 and Res_24_FC > 1 and Sen_24_FC < 1 and Sen_6_FC < 1:
        #         DE_gene_list.append(gene)

        # creat_df = DE_expr_df[DE_expr_df['gene_id'].isin(DE_gene_list) & (DE_expr_df['Ht'] == Ht)]
        # creat_df.to_csv(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_RvSdif_gene.csv',index=False)

if __name__ == "__main__":
    treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat.txt'
    DE_gene_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat_dif_gene.txt'
    R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_fpkm_matrix.csv'
    S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_fpkm_matrix.csv'

    main(treat_file, DE_gene_file, R_count_file, S_count_file)