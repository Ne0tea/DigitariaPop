'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2025-04-19 20:12:52
LastEditors: Ne0tea
LastEditTime: 2025-05-22 20:49:39
'''
import sys
import os
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

pd.set_option('display.max_rows', None)

def process_column_batch(df_block, material_GR50_dic):
    results = []
    for col in df_block.columns:
        col_data = df_block[col]
        counts = col_data.value_counts().reindex([0, 1, 2], fill_value=0)
        mask = (col_data == 2)
        rows_with_2 = col_data.index[mask]
        dict_values = [material_GR50_dic[row] for row in rows_with_2 if row in material_GR50_dic]
        mean_value = np.median(dict_values) if dict_values else 'Nan'
        results.append((col, counts.to_list(), mean_value))
    return results

def main(file_prefix, material_file, max_threads=None, chunk_size=500, outfile='result.txt'):
    file_012 = file_prefix + '.012'
    indv_file = file_012 + '.indv'
    pos_file = file_012 + '.pos'

    # Load index
    with open(indv_file) as f:
        indv_list = [line.strip() for line in f]

    with open(pos_file) as f:
        pos_list = ['_'.join(line.strip().split()) for line in f]

    # Load material data
    material_GR50_dic = {}
    material_list = []
    with open(material_file) as f:
        for line in f:
            key, value = line.strip().split()
            material_GR50_dic[key] = float(value)
            material_list.append(key)

    default_max = min(32, (os.cpu_count() or 1) + 4)
    num_threads = max_threads if max_threads else default_max
    print(f"🧵 Multithread active: {num_threads} threads")
    print(f"📦 Chunk size: {chunk_size} columns")

    # Pre-read only index col
    total_columns = len(pos_list)

    output = open(outfile, 'w')
    for i in range(0, total_columns, chunk_size):
        chunk_cols = list(range(i + 1, min(i + 1 + chunk_size, total_columns + 1)))  # skip col 0
        df_block = pd.read_table(file_012, usecols=[0] + chunk_cols, header=None, index_col=0)
        df_block.index = indv_list
        df_block.columns = pos_list[i:i + len(chunk_cols)]

        hap_GR50_df = df_block.loc[material_list]

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = [executor.submit(process_column_batch, hap_GR50_df, material_GR50_dic)]
            for future in as_completed(futures):
                for col, counts, mean_value in future.result():
                    cur_line = f'{col}' + '\t'.join(map(str, counts)) + f'{mean_value}\n'
                    output.write(cur_line)
    output.close()

if __name__ == "__main__":
    # Usage: python script.py <gene_prefix> <material_file> [threads] [chunk_size]
    file_prefix = sys.argv[1]
    C5_GR50_infv_file = sys.argv[2]
    max_threads = int(sys.argv[3]) if len(sys.argv) > 3 else None
    chunk_size = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    outfile = sys.argv[5]
    main(file_prefix, C5_GR50_infv_file, max_threads, chunk_size, outfile)
