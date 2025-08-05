#statannotayion was run in 3
'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-10-03 20:41:09
LastEditors: Ne0tea
LastEditTime: 2025-07-31 23:27:15
'''
import pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
# template_color={"BAM":"#b71515","ORY":"#e97a0c","POO":"#ffde0a","PAN":"#023e7d","CHL":"#a3cef1"}

def get_diff_expr_gene(expr_count_df, treat_name):
    DE_gene_df = pd.DataFrame(columns=['gene_id','Proj','baseMean','log2FoldChange','lfcSE','pvalue','padj'])
    cur_condition_df = expr_count_df.index.to_series().apply(
        lambda x: 'A' if x.startswith('CK') else ('B' if x.startswith('treatment') else None)
    ).to_frame()
    cur_condition_df.columns = ['condition']

    dds = DeseqDataSet(counts=expr_count_df, metadata=cur_condition_df, fit_type='parametric', design_factors="condition", quiet=False)
    dds.deseq2()
    res = DeseqStats(dds, quiet=True)
    res.summary()
    res_df = res.results_df
    de_res_df = res_df[(res_df['pvalue'] < 0.05)]

    for gene,row in de_res_df.iterrows():
        cur_row = {'gene_id':gene,'Proj':treat_name,'baseMean':row['baseMean'],
                   'log2FoldChange':row['log2FoldChange'],'lfcSE':row['lfcSE'],
                   'pvalue':row['pvalue'],'padj':row['padj']}
        DE_gene_df = DE_gene_df._append(cur_row, ignore_index=True)
    DE_gene_df = DE_gene_df.sort_values(by='gene_id').reset_index(drop=True)

    return DE_gene_df

def read_stringtie_out(file):
    result_df = pd.read_csv(file, header=0)
    result_df = result_df.set_index('gene_id')
    return result_df

def make_treat_dic(treat_file):
    treat_name_SRR={}
    CK_dic = {}
    with open(treat_file, 'r') as gf:
        for i in gf:
            line=i.strip().split()
            srr_id = line[0]
            library_type = line[1].split('|')[1]
            treat_name = line[1].split('|')[0]
            # if srr_id == '15-2-0-3': continue
            if library_type == 'CK': 
                if treat_name.split('_')[0] in CK_dic:
                    CK_dic[treat_name.split('_')[0]].append(srr_id)
                else:
                    CK_dic[treat_name.split('_')[0]] = [srr_id]
                continue
            if treat_name in treat_name_SRR:
                treat_name_SRR[treat_name][library_type].append(srr_id)
            else:
                treat_name_SRR[treat_name]={'CK':[srr_id],'treatment':[]} if library_type == 'CK' else {'CK':[],'treatment':[srr_id]}
    for i in treat_name_SRR:
        cur_type = i.split('_')[0]
        treat_name_SRR[i]['CK'] = CK_dic[cur_type]

    return treat_name_SRR

def make_RvS_treat_dic(treat_file):
    treat_name_SRR={}
    CK_dic = {}
    with open(treat_file, 'r') as gf:
        for i in gf:
            line=i.strip().split()
            srr_id = line[0]
            library_type = line[1].split('|')[1]
            treat_name = line[1].split('|')[0]
            # if srr_id == '15-2-0-3': continue
            if library_type == 'CK': 
                if treat_name.split('_',1)[1] in CK_dic:
                    CK_dic[treat_name.split('_',1)[1]].append(srr_id)
                else:
                    CK_dic[treat_name.split('_',1)[1]] = [srr_id]
                continue
            if treat_name in treat_name_SRR:
                treat_name_SRR[treat_name][library_type].append(srr_id)
            else:
                treat_name_SRR[treat_name]={'CK':[srr_id],'treatment':[]} if library_type == 'CK' else {'CK':[],'treatment':[srr_id]}

    for i in treat_name_SRR:
        cur_type = i.split('_',1)[1]
        treat_name_SRR[i]['CK'] = CK_dic[cur_type]

    return treat_name_SRR

def main(treat_file, R_count_file, S_count_file):

    R_count_df = read_stringtie_out(R_count_file)
    S_count_df = read_stringtie_out(S_count_file)
    all_count_df = pd.concat([R_count_df, S_count_df], axis=1)
    treat_name_SRR = make_RvS_treat_dic(treat_file)
    # treat_name_SRR = make_treat_dic(treat_file)

    DE_gene_df = pd.DataFrame(columns=['gene_id','Proj','baseMean','log2FoldChange','lfcSE','pvalue','padj'])
    for i in treat_name_SRR:
        # if '0h' not in i: continue
        temp_CK_dict = {x:'CK_'+str(idx) for idx, x in enumerate(treat_name_SRR[i]['CK'])}
        temp_Treat_dict = {x:'treatment_'+str(idx) for idx, x in enumerate(treat_name_SRR[i]['treatment'])}
        cur_col_name = temp_CK_dict | temp_Treat_dict

        cur_expr_count_df = all_count_df[list(cur_col_name.keys())]
        cur_expr_count_df = cur_expr_count_df.rename(columns=cur_col_name)
        cur_expr_count_df = cur_expr_count_df.T
        cur_dif_gene_df = get_diff_expr_gene(cur_expr_count_df, i)
        DE_gene_df = pd.concat([DE_gene_df,cur_dif_gene_df],axis=0)
        # print(DE_gene_df[DE_gene_df['gene_id'].isin(['Chr06.4574', 'Chr06.2856'])])
    # print(DE_gene_df)
    # DE_gene_df.to_csv(treat_file.replace('.txt','_dif_gene.txt'),sep='\t',index=False)
    DE_gene_df.to_csv(treat_file.replace('.txt','_RvS_DEGs.csv'),index=False)

if __name__ == "__main__":
    # get 0h control
    # treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat.txt'
    # get S control
    treat_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\Dsan_treat_RvS.txt'
    R_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\21-17-resistance_gene_count_matrix.csv'
    S_count_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_herbicide_DE\15-2-sensitive_gene_count_matrix.csv'

    main(treat_file, R_count_file, S_count_file)