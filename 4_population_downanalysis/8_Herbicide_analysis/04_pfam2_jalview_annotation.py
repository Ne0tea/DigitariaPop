'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-03-21 10:47:08
LastEditors: Ne0tea
LastEditTime: 2025-03-21 15:02:09
'''

pfam_file = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\Herbicide_analysis\Dsan_ALS_pfam.anno'

start_site_list, end_site_list, pfam_list = [], [], []
with open(pfam_file,'r') as pfamf:
    for i in pfamf:
        pfam_name = i.split('\t')[5]
        start_site, end_site = int(i.split('\t')[6]), int(i.split('\t')[7])
        pfam_list.append(pfam_name)
        start_site_list.append(start_site)
        end_site_list.append(end_site)

color_list = ['83c5be', '006d77', 'ffddd2', 'e29578', '0109119']
with open('Jalview_annotation.txt','w') as jalann:
    jalann.write('JALVIEW_ANNOTATION\n')
    jalann.write('SEQUENCE_REF\tChr16.1603.mRNA1\n')
    jalann.write('NO_GRAPH\tPFAM\tPFAM\t')
    switch = False
    for i in range(1,1992):
        if i in start_site_list:
            switch = True
            pfam_name = pfam_list[start_site_list.index(i)]
            cur_color = color_list[pfam_list.index(pfam_name)]
            jalann.write(f'E,{pfam_name},[{cur_color}]|')
            continue
        if i in end_site_list:
            switch = False
            jalann.write(f'E,{pfam_name},[{cur_color}]|')
            continue
        if switch:
            jalann.write(f'E,[{cur_color}]|')
        else:
            jalann.write(f'|')
        # jalann.write('%s\t%s\t%s\t%s\n'%(start_site_list[i], end_site_list[i], pfam_list[i], pfam_list[i]))
    jalann.write('\nCOLOUR\tPFAM\t000000\n')