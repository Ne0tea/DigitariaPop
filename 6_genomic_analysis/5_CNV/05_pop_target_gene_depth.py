'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-15 23:59:34
LastEditors: Ne0tea
LastEditTime: 2025-03-17 13:44:35
'''
import sys
import pandas as pd
import numpy as np
def read_depth_in_file(depth_file):
    gene_dic = {}
    with open(depth_file, 'r') as dfile:
        for i in dfile:
            line = i.strip().split()
            cur_depth = float(line[7])
            for gene in line[6].split(';'):
                if gene in gene_dic:
                    gene_dic[gene].append(cur_depth)
                else:
                    gene_dic[gene] = [cur_depth]
            material_id = line[5]
    gene_dic = {k:sum(v)/len(v) for k,v in gene_dic.items()}
    return gene_dic, material_id

def main(SG_depth_file_list, outfile_name):
    population_depth_df = pd.DataFrame()
    for i in SG_depth_file_list:
        gene_depth_dic, material_id = read_depth_in_file(i)
        gene_depth_df = pd.DataFrame.from_dict(gene_depth_dic, orient="index", columns=[material_id])
        population_depth_df = pd.concat([population_depth_df, gene_depth_df], axis=1)
    population_depth_df.fillna(0, inplace=True)
    population_depth_df.to_csv(outfile_name, sep="\t")

if __name__ == "__main__":
    files = sys.argv[1:-1]
    outfile = sys.argv[-1]
    main(files, outfile)