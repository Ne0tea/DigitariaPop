import sys
homo_site_file=sys.argv[1]
sign7_snp_loc_file=sys.argv[2]

homo_site_list = []
with open(homo_site_file, 'r') as hf:
    for i in hf:
        line=i.strip().split()
        cur_loc = '_'.join(line)
        homo_site_list.append(cur_loc)

with open(sign7_snp_loc_file, 'r') as sf:
    for i in sf:
        line=i.strip().split()
        cur_loc = '_'.join(line)
        if cur_loc not in homo_site_list:
            print(i.strip())
