'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-20 14:05:08
LastEditors: Ne0tea
LastEditTime: 2025-03-20 16:23:39
'''
import sys
import os

visualize_script = r'/public4/home/huangyj/Digitaria/population_downanalysis/Herbicide_analysis/draw_gene_haplotype.R'
herbicide_file = r'/public4/home/huangyj/Digitaria/population_downanalysis/Herbicide_analysis/Vcf_sample_exam_herbicide.txt'
clade_file = r'/public4/home/huangyj/Digitaria/population_downanalysis/Herbicide_analysis/Vcf_sample_exam_clade.txt'

def main(gene_name, herbicide_gz_vcf, gene_ord_file):
    gene_loc_dic = {}
    with open(gene_ord_file, 'r') as gene_ord_f:
        for i in gene_ord_f:
            line=i.strip().split()
            gene_loc_dic[line[3]] = (line[0], int(line[1]), int(line[2]))

    cur_chr, cur_start, cur_end = gene_loc_dic[gene_name]
    os.makedirs('./{}_haplotype'.format(gene_name), exist_ok=True)
    os.chdir('./{}_haplotype'.format(gene_name))
    os.system('bcftools filter ---regions {}:{}-{} {} > {}.vcf'.format(cur_chr, cur_start, cur_end, herbicide_gz_vcf, gene_name))
    os.system('vcftools --vcf {}.vcf --012 --prefix {}'.format(gene_name, gene_name))

    os.system('source activate /public/home/huangyj/anaconda3/envs/population')
    visuliaze_command = f'Rscript {visualize_script} {herbicide_file} {gene_name}.012.pos {clade_file} {gene_name}.012'
    os.system(visuliaze_command)
    os.system('conda deactivate')

    os.chdir('../')

if __name__ == "__main__":
    gene_name = sys.argv[1]
    herbicide_gz_vcf = sys.argv[2]
    gene_ord_file = sys.argv[3]

    main(gene_name, herbicide_gz_vcf, gene_ord_file)
