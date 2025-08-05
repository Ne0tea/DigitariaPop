'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-07 21:03:32
LastEditors: Ne0tea
LastEditTime: 2025-04-13 13:59:00
'''

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['font.family'] = 'Arial'

color_list = ['#d4a373', '#ccd5ae']
sub_color_set={'SubC':'#003049','SubD':'#d62828','SubE':'#f77f00'}
sub_color_set=['#006060', '#B8DADA','#FFBF3D']

def draw_manhadun(signal_df, outdir):
    # sns.set_theme(style="whitegrid")
    certera_99 = np.percentile(signal_df['ZFst'], 99)
    certera_95 = np.percentile(signal_df['ZFst'], 95)
    plot = sns.relplot(data=signal_df, x='SNP_idx', y='ZFst', hue='CHROM', 
                       height=4, aspect=4, 
                       palette = sub_color_set,linewidth=0,s=5, legend=None)
    chrom_df = signal_df.groupby('CHROM')['SNP_idx'].median()
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)

    # 添加各种标注，作为阈值的判断线，调整图像Y轴范围
    plot.ax.set_xlabel('CHROM')
    plot.ax.axhline(y=certera_95, linewidth=2, linestyle="--", color="#003049")
    plot.ax.axhline(y=certera_99, linewidth=2, linestyle="--", color="#780000")
    plot.ax.set_ylim(-1, 10)

    plt.subplots_adjust(top=0.95,bottom=0.1,
                        left = 0.1, right = 0.95)
    # plt.savefig('manhattan.png')
    # plt.show()
    plt.savefig(os.path.join(outdir, r'RvS_Fst.pdf'), format='pdf')
    top_95_percent = signal_df[signal_df['ZFst'] > certera_95]
    return top_95_percent

def draw_pi(signal_df, outdir):
    # sns.set_theme(style="whitegrid")
    certera_99 = np.percentile(signal_df['W'], 99)
    certera_95 = np.percentile(signal_df['W'], 95)
    plot = sns.relplot(data=signal_df, x='region_start', y='W', hue='CHROM',
                       height=4, aspect=4, linewidth=0,s=5,
                       palette = sub_color_set, legend=None)
    chrom_df = signal_df.groupby('CHROM')['region_start'].median()
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)

    plot.ax.set_xlabel('')
    plot.ax.set_ylabel('PiS/PiR')
    plot.ax.axhline(y=certera_95, linewidth=2, linestyle="--", color="#003049")
    plot.ax.axhline(y=certera_99, linewidth=2, linestyle="--", color="#780000")
    # plot.ax.set_ylim(0, 10)

    plt.subplots_adjust(top=0.95,bottom=0.1,
                        left = 0.1, right = 0.95)
    plt.savefig(os.path.join(outdir, r'RvS_Pi.pdf'), format='pdf')
    # plt.show()
    top_95_percent = signal_df[signal_df['W'] > certera_95]
    return top_95_percent

def draw_manhadun_dxy(signal_df, outdir):
    # sns.set_theme(style="whitegrid")
    certera_99 = np.percentile(signal_df['dxy_Resistant_Sensitive'], 99)
    # print(signal_df['dxy_Resistant_Sensitive'], certera_99)
    certera_95 = np.percentile(signal_df['dxy_Resistant_Sensitive'], 95)
    plot = sns.relplot(data=signal_df, x='SNP_idx', y='dxy_Resistant_Sensitive', hue='scaffold', 
                       height=4, aspect=4, 
                       palette = sub_color_set,linewidth=0,s=5, legend=None)
    chrom_df = signal_df.groupby('scaffold')['SNP_idx'].median()

    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)

    plot.ax.set_xlabel('scaffold')
    plot.ax.axhline(y=certera_95, linewidth=2, linestyle="--", color="#003049")
    plot.ax.axhline(y=certera_99, linewidth=2, linestyle="--", color="#780000")
    # plot.ax.set_ylim(-1, 10)

    plt.subplots_adjust(top=0.95,bottom=0.1,
                        left = 0.1, right = 0.95)
    # plt.savefig('manhattan.png')
    # plt.show()
    plt.savefig(outdir + r'RvS_dxy.pdf', format='pdf')
    top_95_percent = signal_df[signal_df['dxy_Resistant_Sensitive'] > certera_95]
    return top_95_percent

def main(signal_file, R_file, S_file, dxy_signal_file, GWAS_candidate_loc_file):
    with open(GWAS_candidate_loc_file, 'r') as herbf:
        herb_loc_dic = {}
        for i in herbf:
            line = i.strip().split('\t')
            if int(line[0]) < 10:
                chr_name = 'Chr0' + line[0]
            else:
                chr_name = 'Chr' + line[0]
            if chr_name not in herb_loc_dic:
                herb_loc_dic[chr_name] = [int(line[2])]
            else:
                herb_loc_dic[chr_name].append(int(line[2]))

    signal_df = pd.read_table(signal_file,header=0)

    mean = signal_df['WEIGHTED_FST'].mean()
    std = signal_df['WEIGHTED_FST'].std()
    signal_df['ZFst'] = (signal_df['WEIGHTED_FST'] - mean) / std
    signal_df = signal_df[signal_df['CHROM'].str.contains('Chr')]
    signal_df["SNP_idx"]=signal_df.index
    out_dir = os.path.dirname(signal_file)
    top_region_FST = draw_manhadun(signal_df, out_dir)
    top_region_FST.to_csv(signal_file + r'.top_region_Fst.txt',sep='\t',index=False)

    R_df = pd.read_table(R_file,header=0)
    R_df = R_df[R_df['CHROM'].str.contains('Chr')]
    R_df.columns = ['CHROM','START', 'END', 'R_VAR','R_pi']
    S_df = pd.read_table(S_file,header=0)
    S_df = S_df[S_df['CHROM'].str.contains('Chr')]
    S_df.columns = ['CHROM','START', 'END', 'S_VAR','S_pi']
    All_pi_df = pd.merge(R_df,S_df,on=['CHROM','START', 'END'])
    All_pi_df['W'] = All_pi_df['S_pi'] / All_pi_df['R_pi']
    All_pi_df["region_start"]=All_pi_df.index
    out_dir = os.path.dirname(R_file)
    top_region_Pi = draw_pi(All_pi_df, out_dir)
    top_region_Pi.to_csv(R_file.replace('R','RvS') + r'.top_region_Pi.txt',sep='\t',index=False)

    with open(R_file.replace('R','RvS') + r'.top_region_Pi.txt', 'r') as trpf:
        for i in trpf:
            line = i.strip().split('\t')
            if line[0] not in herb_loc_dic: continue
            for snp_loc in herb_loc_dic[line[0]]:
                if int(line[1]) <= int(snp_loc) and int(line[2]) >= int(snp_loc):
                    print(snp_loc, i.strip())

    with open(signal_file + r'.top_region_Fst.txt', 'r') as trpf:
        for i in trpf:
            line = i.strip().split('\t')
            if line[0] not in herb_loc_dic: continue
            for snp_loc in herb_loc_dic[line[0]]:
                if int(line[1]) <= int(snp_loc) and int(line[2]) >= int(snp_loc):
                    print(snp_loc, i.strip())

    # dxy_df = pd.read_csv(dxy_signal_file,header=0)
    # dxy_df["SNP_idx"]=dxy_df.index
    # dxy_df.dropna(inplace=True)
    # out_dir = os.path.dirname(dxy_signal_file)
    # top_region_dxy = draw_manhadun_dxy(dxy_df, out_dir)
    # top_region_dxy.to_csv(dxy_signal_file.replace('.Fst.Dxy.pi','RvS') + r'.top_region_Dxy.txt',sep='\t',index=False)

if __name__ == "__main__":
    Fst_signal_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Dsan_RvS_NE.100k50k.windowed.weir.fst'
    dxy_signal_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.herbicide.50k.Fst.Dxy.pi.csv'
    pi_R_signal_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Dsan_R_NE.100k50k.windowed.pi'
    pi_S_signal_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Selective_sweep\Dsan_S_NE.100k50k.windowed.pi'

    GWAS_candidate_loc_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_all.merge.gatk3.fd.emmax.3.covkin-TopTenPC.assoc.txt.sign7.txt'
    main(Fst_signal_file, pi_R_signal_file, pi_S_signal_file, dxy_signal_file, GWAS_candidate_loc_file)