'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-19 14:00:45
LastEditors: Ne0tea
LastEditTime: 2025-04-11 20:32:29
'''
import pandas as pd
import os
import sys

def find_overlapping_genes(row, snp_dict):
    chr = row['Chr']
    start = int(row['START'])
    end = int(row['END'])
    overlapping_genes = []
    if chr in snp_dict:
        for loc, snp_info in snp_dict[chr].items():
            snp_loc = int(loc)
            if (start <= snp_loc) and (end >= snp_loc):
                overlapping_genes.append(str(snp_loc))
    return ', '.join(overlapping_genes) if overlapping_genes else ''

def main(snp_file, GWAS_type, snp_eff_file):
    filename = os.path.basename(snp_file)
    snp_eff_dic = {}
    with open(snp_eff_file, 'r') as gene_ord_f:
        for i in gene_ord_f:
            if i.startswith('#'): continue
            line=i.strip().split()
            if line[0] in snp_eff_dic:
                snp_eff_dic[line[0]][line[1]]= line[2:]
            else:
                snp_eff_dic[line[0]] = {line[1]: line[2:]}

    if GWAS_type == 'emmax':
        snp_df = pd.read_table(snp_file, sep = ' ',  names=['CHR', 'CHR_BP_REF_ALT', 'BP', 'P', 'BP_ID'])
    else:
        snp_df = pd.read_table(snp_file, sep = '\t',  names=['CHR', 'CHR_BP_REF_ALT', 'BP', 'Nmiss', 'Refalleo','Altalleo','Af','Beta','Se','Logl_H1','L_remle','P'])
    snp_df['START'] = snp_df['BP'] - 130000
    snp_df['END'] = snp_df['BP'] + 130000
    snp_df['Chr'] = snp_df['CHR'].apply(lambda x: 'Chr' + str(x) if int(x) >= 10 else 'Chr0'+str(x))
    snp_df['SNP'] = snp_df.apply(lambda row: find_overlapping_genes(row, snp_eff_dic), axis=1)

    with open(snp_file + r'_candidate_SNPs.txt.txt', 'w') as of:
        for row, rowdf in snp_df.iterrows():
            if rowdf['SNP']:
                for snp in rowdf['SNP'].split(', '):
                    of.write('\t'.join([rowdf['Chr'], rowdf['CHR_BP_REF_ALT'], snp, '\t'.join(snp_eff_dic[rowdf['Chr']][snp]) ])+'\n')
            else:
                of.write('\t'.join([rowdf['Chr'], rowdf['CHR_BP_REF_ALT'], '', '\t'.join(snp_eff_dic[rowdf['Chr']][snp])])+'\n')

if __name__ == "__main__":
    snp_file = sys.argv[1]
    GWAS_type = sys.argv[2]

    snp_eff_file = r'/public4/home/huangyj/Digitaria/population_downanalysis/vcf_filt/all_gatk3/Dsan_all.merge.gatk3.fd.nohomo.recode_fded_SNP.eff.annotations.txt'
    main(snp_file, GWAS_type, snp_eff_file)
