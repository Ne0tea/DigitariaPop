'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-12 01:48:24
LastEditors: Ne0tea
LastEditTime: 2025-04-13 16:32:37
'''
import sys

def main(snp_var_file, target_gene_file):
    gene_name_list=[]
    with open(target_gene_file, 'r') as target_gene_f:
        for i in target_gene_f:
            gene_name_list.append(i.strip())

    with open(snp_var_file, 'r') as snp_var_f:
        for i in snp_var_f:
            if i.startswith('#'): continue
            line=i.strip().split()
            fields = line[7].split(',')
            for field in fields:
                field = field.split('|')
                var_type = field[1]
                var_serve = field[2]
                var_gene = field[3]
                if var_gene not in gene_name_list:
                    continue
                if 'upstream' in var_type or 'downstream' in var_type:
                    gene_dis = int(field[14])
                    if gene_dis < 5000:
                        print(i.strip())
                        break
                if 'missense' in var_type or 'synonymous' in var_type or 'intron_variant' in var_type:
                    print(i.strip())
                    break
                if var_serve == 'MODERATE' or var_serve == 'HIGH':
                    print(i.strip())
                    break
if __name__ == "__main__":
    snp_eff_file = sys.argv[1]
    target_gene_file = sys.argv[2]

    # snp_eff_file = r'F:\GWAS\test.snp'
    main(snp_eff_file, target_gene_file)
