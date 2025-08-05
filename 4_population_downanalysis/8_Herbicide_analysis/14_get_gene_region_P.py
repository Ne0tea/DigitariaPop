'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-22 21:11:50
LastEditors: Ne0tea
LastEditTime: 2025-05-22 22:20:43
'''
import sys
def main(gene_bed, log_file, gene_list_file):
    signal_dic = {}
    with open(log_file, 'r') as lgf:
        for i in lgf:
            if i.endswith('p_wald\n'):
                continue
            line=i.strip().split()
            cur_chr = 'Chr' + str(line[0]) if int(line[0])>=10 else 'Chr0' + str(line[0])
            cur_loc = int(line[2])
            if cur_chr in signal_dic:
                signal_dic[cur_chr][cur_loc] = line[-1]
            else:
                signal_dic[cur_chr] = {cur_loc:line[-1]}

    gene_dic = {}
    with open(gene_bed, 'r') as gbf:
        for i in gbf:
            line = i.strip().split()
            gene_dic[line[3]] = (line[0], int(line[1]), int(line[2]))

    with open(gene_list_file, 'r') as of:
        for i in of:
            line = i.strip().split()
            cur_chr, cur_s, cur_e = gene_dic[line[0]]
            for loc in signal_dic[cur_chr]:
                if loc >= cur_s and loc <= cur_e:
                    print(line[0], cur_s, cur_e, loc, signal_dic[cur_chr][loc], sep='\t')

if __name__ == "__main__":
    gene_bed = sys.argv[1]
    log_file = sys.argv[2]
    gene_list_file = sys.argv[3]
    main(gene_bed, log_file, gene_list_file)
