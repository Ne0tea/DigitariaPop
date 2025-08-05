import pandas as pd
import argparse

def get_consensus(infile, outfile,mnumber):
    df_Summary = pd.read_table(infile, sep=",") #name start end ave.score most.freq.value.N consensus.primary consensus.count width fasta.name class consensus.secondary repeats.identified
    consensus_repeat_max_num_list = df_Summary[(df_Summary['most.freq.value.N'] > 150) & (df_Summary['most.freq.value.N'] < 650)]['repeats.identified'].nlargest(int(mnumber))
    #consensus_repeat_max_num_list = df_Summary['repeats.identified'].nlargest(int(mnumber))
    o = open(outfile, 'w') 
    o.write("name,length,seq\n")
    for consensus_repeat_max_num in consensus_repeat_max_num_list:
        repeat_length = int(df_Summary[df_Summary['repeats.identified'] == consensus_repeat_max_num]['most.freq.value.N'].iloc[0])
        seq = df_Summary.loc[df_Summary['repeats.identified'] == consensus_repeat_max_num]['consensus.secondary'].iloc[0]
        o.write("CEN" + str(repeat_length) + '_'+str(consensus_repeat_max_num) +"," + str(repeat_length) + "," + str(seq)+'\n')
    o.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='summary_repeats', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output', dest='output', required=True)
    parser.add_argument('-m', '--maxnumber', help='output monomer number', dest='mnumber', required=True)
    args = parser.parse_args()
    get_consensus(args.input, args.output,args.mnumber)

if __name__ == '__main__':
    main()
