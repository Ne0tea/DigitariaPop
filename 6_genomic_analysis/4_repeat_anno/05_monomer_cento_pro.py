'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-07-25 16:07:19
LastEditors: Ne0tea
LastEditTime: 2024-07-25 20:49:44
'''
import pandas as pd
import argparse
def get_monomer_centromere(CR_region_file, repeat_file, prefix, outfile, CR_length):
    df_cent = pd.read_table(CR_region_file, sep='\t', low_memory=False,names=['Material','chr','horclass','start','end','length'])
    df_all = pd.read_table(repeat_file, sep=',', low_memory=False) # start,end,width,seq,strand,class,region.name,seq.name,edit.distance,repetitiveness
    df_all = df_all.map(str)
    df_all = df_all.map(remove_quotes)
    df_all = df_all.map(custom_conversion)
    chr_list = df_all['seq.name'].unique()
    CR_list = df_cent['chr'].unique()
    df_all['start'] = df_all['start'].astype(float)
    df_all['end'] = df_all['end'].astype(float)
    o = open(outfile, 'w')
    #Dsan    Chr01   CEN159  9378001 9568699 190699  0.9793181925442713
    o.write("Material" + "\t" + "Chr" + "\t" + "Horclass" + "\t" + "Start" + "\t" + "End" + "\t" + "Length" + "\t" + "Rproportion(%)" + "\n")

    for each_chr in chr_list:
        if each_chr not in CR_list:
            o.write(str(prefix) + "\t" + str(each_chr) + "\t" + "no centromere" + "\n")
        else:
            df_length_centromere = df_cent[(df_cent['Material'] == str(prefix)) & (df_cent['chr'] == str(each_chr)) & (df_cent['length'] > CR_length)].reset_index()
            for index, row in df_length_centromere.iterrows():
                cen_s = row['start']
                cen_e = row['end']
                length_centromere =  int(cen_e) - int(cen_s) + 1
                #df_CEN = df_all[(df_all['seq.name'] == str(each_chr)) & (df_all['class'].str.contains('CEN')) & \
                df_CEN = df_all[(df_all['seq.name'] == str(each_chr)) & (df_all['width'] > 150) & \
                                (df_all['start'] >= cen_s) & (df_all['end'] <= cen_e)].reset_index()
                length_total=df_CEN['width'].sum()
                pro = length_total / length_centromere
                o.write(str(prefix) + "\t" + str(each_chr)+ "\t" + str(row['horclass'])+ "\t" \
                            + str(cen_s)+ "\t" + str(cen_e)+ "\t" + str(length_centromere) + "\t" + str(pro) + "\n")
    o.close()

def remove_quotes(text):
    return text.replace('"', '')

def custom_conversion(cell_value):
    try:
        return int(cell_value)  # 尝试将字符串转换为整数
    except ValueError:
        return cell_value  # 如果无法转换，保留原字符串

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='prefix', dest='pre', required=True)
    parser.add_argument('-cr', '--cregion', help='CR region file', dest='cregion', required=True)
    parser.add_argument('-mr', '--mrepeat', help='monomer repeat', dest='mrepeat', required=True)
    parser.add_argument('-o', '--output', help='output', dest='outfile', required=True)
    parser.add_argument('-crl', '--crlength', help='CR length', dest='CRlength', default=5000)

    args = parser.parse_args()
    get_monomer_centromere(args.cregion, args.mrepeat, args.pre, args.outfile,args.CRlength)

if __name__ == '__main__':
    main()

