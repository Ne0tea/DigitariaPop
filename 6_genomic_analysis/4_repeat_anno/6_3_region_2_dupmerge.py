import sys

an_fai = sys.argv[1]
bed_file = sys.argv[2]

chr_length_dic={}
with open(an_fai,'r') as faif:
    for i in faif:
        line=i.strip().split()
        chr_length_dic[line[0]]=int(line[1])

outfile=open(bed_file + '.modi', 'w')
with open(bed_file, 'r') as ttf:
    for i in ttf:
        line=i.strip().split()
        cur_chr=line[0]
        if cur_chr not in chr_length_dic:
            continue
        cur_chr_length=chr_length_dic[cur_chr]
        cur_line='\t'.join([cur_chr, str(cur_chr_length), line[1], line[2], str(int(line[2])-int(line[1])), line[3].split('-')[1]])
        outfile.write(cur_line+'\n')

outfile.close()
