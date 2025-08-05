import argparse
import pandas as pd


def get_centromere_repeat(allrepeat_file, prefix, outfile, region_min):
    df = pd.read_table(allrepeat_file, sep=',') # start,end,width,seq,strand,class,region.name,seq.name,edit.distance,repetitiveness
    chr_list = df.drop_duplicates(subset=['seq.name'], keep='last', ignore_index=True)['seq.name'].tolist()
    CEN_type = df.drop_duplicates(subset=['class'], keep='last', ignore_index=True)['class'].tolist()

    ### get rid of nan value
    CEN_type_new = [elem if not pd.isnull(elem) else None for elem in CEN_type]
    while None in CEN_type_new:
        CEN_type_new.remove(None)

    o = open(outfile, 'w')
    # o.write("Accession" + "\t" + "Chr" + "\t" + "CEN_type" + "\t" + "Start" + "\t" + "End" + "\t" + "Length" + "\n")
    for i in chr_list:
        for j in CEN_type_new:
            df_chr = df[(df['seq.name'] == str(i)) & (df['class'] == str(j)) & (df['width'] >= 140) & (df['width'] <= 170)].reset_index(drop=True)
            if df_chr.empty:
                continue
            circle = 0
            line_num = 0
            next_line = 1
            while next_line < df_chr.shape[0]:
                if (int(sum((df_chr.iloc[line_num:(next_line + 1), 2])))) <= (int(df_chr.iloc[next_line, 1]) - int(df_chr.iloc[line_num, 0]) + 1) <= (int(sum((df_chr.iloc[line_num:(next_line + 1), 2]))) + 150000) :
                    circle += 1
                    line_num += 1
                    next_line += 1
                else:
                    start_line0 = line_num - circle
                    end_line0 = next_line - 1
                    if end_line0 - start_line0 > 1:
                        start0 = df_chr.loc[start_line0, "start"]
                        end0 = df_chr.loc[end_line0, "end"]
                        length0 = end0 -start0 + 1
                        o.write(str(prefix) + "\t" + str(i) + "\t" + str(j) + "\t" + str(int(start0)) + "\t" + str(int(end0)) + "\t" + str(int(length0)) + "\n")
                        circle = 0
                        line_num += 1
                        next_line += 1
                    else:
                        circle = 0
                        line_num += 1
                        next_line += 1
            start_line = line_num - circle
            end_line = next_line - 1
            start = df_chr.loc[start_line, "start"]
            end = df_chr.loc[end_line, "end"]
            length = end - start + 1
            o.write(str(prefix) + "\t" + str(i) + "\t" + str(j) + "\t" + str(int(start)) + "\t" + str(int(end)) + "\t" + str(int(length)) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='all repeats file from TRASH', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    parser.add_argument('-p', '--prefix', help='prefix for outfile', dest='prefix', required=True)
    parser.add_argument('--min_region_length', help='min length for repeat region', dest='region', default=5000)
    args = parser.parse_args()
    get_centromere_repeat(args.input, args.prefix, args.output, args.region)

if __name__ == '__main__':
    main()
