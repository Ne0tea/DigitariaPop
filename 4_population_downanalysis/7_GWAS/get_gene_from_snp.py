'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-19 14:00:45
LastEditors: Ne0tea
LastEditTime: 2025-03-19 15:02:02
'''
import numpy as np
import pandas as pd
import sys
import os

def find_overlapping_genes(row, gene_dict):
    chr = row['Chr']
    start = int(row['START'])
    end = int(row['END'])
    overlapping_genes = []
    if chr in gene_dict:
        for gene, (g_start, g_end) in gene_dict[chr].items():
            g_start = int(g_start)
            g_end = int(g_end)
            if (start <= g_end) and (end >= g_start):
                overlapping_genes.append(gene)
    return ', '.join(overlapping_genes) if overlapping_genes else ''


def main(snp_file, GWAS_type, gene_ord_file):
    filename = os.path.basename(snp_file)
    gene_loc_dic = {}
    with open(gene_ord_file, 'r') as gene_ord_f:
        for i in gene_ord_f:
            line=i.strip().split()
            if line[0] in gene_loc_dic:
                gene_loc_dic[line[0]][line[3]]= (int(line[1]), int(line[2]))
            else:
                gene_loc_dic[line[0]] = {line[3]: (int(line[1]), int(line[2]))}

    if GWAS_type == 'emmax':
        snp_df = pd.read_table(snp_file, sep = ' ',  names=['CHR', 'CHR_BP_REF_ALT', 'BP', 'P', 'BP_ID'])
    else:
        snp_df = pd.read_table(snp_file, sep = '\t',  names=['CHR', 'CHR_BP_REF_ALT', 'BP', 'Nmiss', 'Refalleo','Altalleo','Af','Beta','Se','Logl_H1','L_remle','P'])
    snp_df['START'] = snp_df['BP'] - 130000
    snp_df['END'] = snp_df['BP'] + 130000
    snp_df['Chr'] = snp_df['CHR'].apply(lambda x: 'Chr' + str(x) if int(x) >= 10 else 'Chr0'+str(x))
    snp_df['Genes'] = snp_df.apply(lambda row: find_overlapping_genes(row, gene_loc_dic), axis=1)

    candidate_genes = set([ y for x in snp_df['Genes'].to_list() for y in x.split(', ') if x])
    with open(snp_file + r'_candidate_genes.txt', 'w') as of:
        for i in sorted(list(candidate_genes)):
            of.write(i+'\n')

if __name__ == "__main__":
    snp_file = sys.argv[1]
    GWAS_type = sys.argv[2]
    # snp_file = r'F:\GWAS\emmax\Dsan_all.merge.gatk3.vcf.filtered.QC.Inbreeding.LD.emmax.GR50.significant5.txt'
    # snp_file = r'F:\GWAS\emmax\Dsan_all.merge.gatk3.vcf.filtered.QC.Inbreeding.LD.emmax.GR50.significant7.txt'
    # snp_file = r'F:\GWAS\emmax\Dsan_all.merge.gatk3.vcf.filtered.QC.Inbreeding.LD.emmax.GR90.significant8.txt'
    # snp_file = r'F:\GWAS\emmax\Dsan_all.merge.gatk3.vcf.filtered.QC.Inbreeding.LD.emmax.GR90.significant10.txt'

    gene_ord_file = r'/public4/home/huangyj/Digitaria/related_ref/Dsan_V3.ord.gene.bed'
    main(snp_file, GWAS_type, gene_ord_file)
