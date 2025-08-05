'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-05-20 14:04:56
LastEditors: Ne0tea
LastEditTime: 2025-07-10 23:44:21
'''
from collections import defaultdict
def calculate_chromosome_averages(file_path):
    chrom_sums = defaultdict(lambda: [0, 0])
    with open(file_path) as f:
        for line in f:
            if line.startswith('#'):  # 跳过注释行
                continue
            parts = line.strip().split()
            if len(parts) >= 5:
                chrom = parts[0]
                value = float(parts[4])
                chrom_sums[chrom][0] += value
                chrom_sums[chrom][1] += 1

    chrom_avgs = {}
    for chrom, (total, count) in chrom_sums.items():
        chrom_avgs[chrom] = total / count if count != 0 else 0
    return chrom_avgs

def compare_files(file1_path, file2_path, subg_chr_file):
    file1_avgs = calculate_chromosome_averages(file1_path)
    file2_avgs = calculate_chromosome_averages(file2_path)

    subg_dic = {}
    with open(subg_chr_file,'r') as subg_file:
        for i in subg_file:
            subg_dic[i.strip().split('\t')[0]] = i.strip().split('\t')[2]

    file2_data = {}
    with open(file2_path) as f2:
        for line in f2:
            if line.startswith('#'):  # 跳过注释行
                continue
            parts = line.strip().split()
            if len(parts) >= 5:  # 确保有至少5列
                key = (parts[0], parts[1], parts[2])  # 前三列作为键
                file2_data[key] = float(parts[4])  # 第五列作为值

    results = []
    with open(file1_path) as f1:
        for line in f1:
            if line.startswith('#'):  # 跳过注释行
                continue
            parts = line.strip().split()
            if parts[0] not in subg_dic:
                continue
            cur_subg = subg_dic[parts[0]]
            cur_subg_avg_in_file1 = file1_avgs[parts[0]]
            cur_subg_avg_in_file2 = file2_avgs[parts[0]]
            if len(parts) >= 5:  # 确保有至少5列
                key = (parts[0], parts[1], parts[2])
                value1 = float(parts[4])
                # print(value1, value2)
                if key in file2_data:
                    value2 = file2_data[key]
                    if cur_subg == 'SubC':
                        if value2 > value1:
                            results.append((line.strip(), value1, value2))
                    elif cur_subg == 'SubD' or cur_subg == 'SubE':
                        if value2 < value1:
                            results.append((line.strip(), value1, value2))
    return results

if __name__ == "__main__":
    #file1 is SubC
    file1=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Phylogeny\YZGJ2.sorted.bam.100k50k.depth'
    #file2 is SubD SubE
    file2=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\Phylogeny\DZ2.sorted.bam.100k50k.depth'
    subg_chr_file=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Geomic_analysis\0_Dsan_chr_trans'

    matched_lines = compare_files(file1, file2, subg_chr_file)
    for region, v1, v2 in matched_lines:
        print(region, v1, v2)
