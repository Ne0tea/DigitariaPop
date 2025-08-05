'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-15 23:59:34
LastEditors: Ne0tea
LastEditTime: 2025-03-16 00:17:57
'''
import sys

def main(SG_depth_file, target_gene_file, out_file):
    SG_depth_listr = []
    with open(SG_depth_file, 'r') as f:
        for line in f:
            if 'readCount' in line: continue
            line = line.strip()
            cur_depth = float(line[4])
            SG_depth_listr.append(cur_depth)
    average_depth = sum(SG_depth_listr)/len(SG_depth_listr)
    output = open(out_file,'w')
    with open(target_gene_file, 'r') as f:
        for line in f:
            if 'readCount' in line: continue
            line = line.strip()
            gene_depth = float(line.split('\t')[4]) / average_depth
            output.write(line+'\t'+str(gene_depth)+'\n')
    output.close()

if __name__ == "__main__":
    SG_depth_file = sys.argv[1]
    target_gene_file = sys.argv[2]
    out_file = sys.argv[3]
    main(SG_depth_file, target_gene_file, out_file)