'''
Descripttion: generate itol annotation file
Author: Ne0tea
version: 
Date: 2023-05-14 21:56:06
LastEditors: Ne0tea
LastEditTime: 2025-04-24 12:04:41
'''
import sys

group_color={'Group_1':'#006d77','Group_2':'#42999b','Group_3':'#83c5be','Group_4':'#b8dedc','Group_5':'#edf6f9',
                'Group_6':'#f6eae6','Group_7':'#ffddd2','Group_8':'#f1b9a5','Group_12':'#e29578'}
type_color={'sympatric':'#f4a261','allopatric':'#e9c46a','C4':'#2a9d8f','Outgroup':'#264653'}

def write_group_itol_anno(pop_file, config_annotation_file):
    anno_out = open(config_annotation_file, 'w')
    anno_out.write("DATASET_COLORSTRIP\nSEPARATOR\tTAB\n")
    anno_out.write("BORDER_WIDTH\t0.5\n")
    anno_out.write("COLOR\t#bebada\n")
    anno_out.write("DATASET_LABEL\tGroup number\n")

    color_line="\t".join(group_color.values())
    label_line="\t".join(group_color.keys())
    shape_line="\t".join(str(1)*len(group_color))

    anno_out.write("LEGEND_COLORS\t"+color_line+"\n")
    anno_out.write("LEGEND_LABELS\t"+label_line+"\n")
    anno_out.write("LEGEND_SHAPES\t"+shape_line+"\n")
    anno_out.write("LEGEND_TITLE\tClade\n")
    anno_out.write("MARGIN\t5\n")
    anno_out.write("STRIP_WIDTH\t25\n")
    anno_out.write("DATA\n")

    with open(pop_file, 'r') as pf:
        for i in pf:
            line=i.strip().split()
            cur_id = line[0]
            cur_clade=line[2]
            cur_color = group_color[cur_clade]
            c_line="\t".join([cur_id,cur_color,cur_clade])+"\n"
            anno_out.write(c_line)
        anno_out.close()

    print("Annotation successfully stored in "+config_annotation_file)

def write_type_itol_anno(pop_file, config_annotation_file):
    anno_out = open(config_annotation_file, 'w')
    anno_out.write("DATASET_COLORSTRIP\nSEPARATOR\tTAB\n")
    anno_out.write("BORDER_WIDTH\t0.5\n")
    anno_out.write("COLOR\t#bebada\n")
    anno_out.write("DATASET_LABEL\tSample Type\n")

    color_line="\t".join(type_color.values())
    label_line="\t".join(type_color.keys())
    shape_line="\t".join(str(1)*len(type_color))

    anno_out.write("LEGEND_COLORS\t"+color_line+"\n")
    anno_out.write("LEGEND_LABELS\t"+label_line+"\n")
    anno_out.write("LEGEND_SHAPES\t"+shape_line+"\n")
    anno_out.write("LEGEND_TITLE\tClade\n")
    anno_out.write("MARGIN\t5\n")
    anno_out.write("STRIP_WIDTH\t25\n")
    anno_out.write("DATA\n")

    used_id = []
    with open(pop_file, 'r') as pf:
        for i in pf:
            line=i.strip().split()
            cur_id = line[0]
            if cur_id in used_id:
                continue
            cur_type=line[1]
            cur_color = type_color[cur_type]
            c_line="\t".join([cur_id,cur_color,cur_type])+"\n"
            anno_out.write(c_line)
        anno_out.close()

    print("Annotation successfully stored in "+config_annotation_file)

def main(pop_file):
    group_annotation_file = pop_file + '_Group_itol_anno.txt'
    type_annotation_file = pop_file + '_Type_itol_anno.txt'
    print(f"Pop file: {pop_file}")
    print(f"Output group config file: {group_annotation_file}")
    print(f"Output type config file: {type_annotation_file}")

    write_group_itol_anno(pop_file, group_annotation_file)
    write_type_itol_anno(pop_file, type_annotation_file)

if __name__ == "__main__":
    pop_file = sys.argv[1]
    main(pop_file)
