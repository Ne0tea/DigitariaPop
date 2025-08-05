'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-15 23:59:34
LastEditors: Ne0tea
LastEditTime: 2025-08-01 23:30:37
'''
import sys
import os
import pandas as pd
import numpy as np
from joypy import joyplot
import matplotlib.pyplot as plt
plt.rcParams[ 'pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams[ 'font.family'] = 'Arial'
ignore_material = ['YZGJ2_1','NT3_1']
joy_color = {'Resistant':'#dda15e','Sensitive':'#ccd5ae'}
Doamin_color = {'ALS':'#ffddd2', 'EPSPS':'#edf6f9', 'UGT':'#e29578', 'AKR':'#f9ebc7', 'ABC':'#ccd5ae', 'NB-ARC':'#d4a373', 'GST':'#83c5be', 'P450':'#006d77'}

sub_color_set={'SubC':'#003049','SubD':'#d62828','SubE':'#f77f00'}
def draw_joyplot(i, depth_df, Vcf_RR_dic, out_dir, prefix):
    # depth_df = depth_df.drop(columns=ignore_material)

    depth_plot_df = pd.melt(
    depth_df,
    id_vars=['geneID'],
    var_name='Vcf_id',
    value_name='Depth'
    )
    depth_plot_df['Type'] = depth_plot_df['Vcf_id'].map(Vcf_RR_dic)
    depth_plot_df.dropna(inplace=True)
    depth_plot_df = depth_plot_df.pivot(
        index=['geneID', 'Vcf_id'],  # 保留不变的列
        columns='Type',
        values='Depth'
    ).reset_index()

    fig,ax = joyplot(depth_plot_df, by='geneID', column=['Resistant', 'Sensitive'], figsize= (4,4), legend=True,
            color=['#dda15e','#ccd5ae'], linewidth=.8, grid=True, alpha=.3, x_range=[-1, 5])

    plt.subplots_adjust(top=0.95,bottom=0.1)
    plt.show()
    # plt.savefig(os.path.join(out_dir,prefix+'_'+str(i*100)+'.pdf'), format='pdf')
    plt.close()

def main(depth_file, herbicide_resistance_file, out_dir, prefix):
    Vcf_RR_dic={}
    C5_material_list = []
    with open(herbicide_resistance_file, 'r') as f:
        for i in f:
            line = i.strip().split()
            Vcf_RR_dic[line[0]] = line[1]
            C5_material_list.append(line[0])

    depth_df = pd.read_table(depth_file)
    depth_df.columns = ["geneID"] + depth_df.columns[1:].tolist()
    depth_df = depth_df[C5_material_list+['geneID']]
    chunks = np.array_split(depth_df, np.ceil(len(depth_df) / 100))

    os.makedirs(out_dir, exist_ok=True)
    for i, chunk in enumerate(chunks):
        print(f"Processing chunk {i+1}:")
        draw_joyplot(i, chunk, Vcf_RR_dic, out_dir, prefix)


if __name__ == "__main__":
    # files = sys.argv[1:-1]
    # outfile = sys.argv[-1]
    # for i in ['UGT','P450','ABC','GST','NB-ARC','ALS','EPSPS']:
    # for i in ['GST','NB-ARC','ALS','EPSPS']:
    for i in ['ALS']:
        depth_file = rf'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\{i}_uniq_ratio.txt'
        herbicide_resistance_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Digitaria_HB_RR.txt'
        out_dir = rf'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\CNV_analysis\{i}_dep_file'
        prefix = i
        main(depth_file, herbicide_resistance_file, out_dir, prefix)