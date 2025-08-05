'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-08-26 13:31:36
LastEditors: Ne0tea
LastEditTime: 2025-04-30 14:29:11
'''
import pandas as pd
from geopy.distance import great_circle
import numpy as np
import os

def within_50km(lat1, lon1, lat2, lon2,distance_threshold_km):
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    return great_circle(point1, point2).kilometers <= distance_threshold_km

def main(C4_df,C5_df,distance_threshold_km,out_dir):
    C4_sympatric = pd.DataFrame(columns=C4_df.columns)
    C5_sympatric = pd.DataFrame(columns=C5_df.columns)

    # 遍历 A 中的每一行，筛选出 B 和 C 中距离50公里以内的点，并合并到新的数据框中
    for index, row in C4_df.iterrows():
        lat, lon = row['Lat'], row['Lon']

        C4_within_50km = C4_df[C4_df.apply(lambda x: within_50km(lat, lon, x['Lat'], x['Lon'],distance_threshold_km), axis=1)]
        C5_within_50km = C5_df[C5_df.apply(lambda x: within_50km(lat, lon, x['Lat'], x['Lon'],distance_threshold_km), axis=1)]
        if len(C4_within_50km) >= 3 and len(C5_within_50km) >= 3:
            print(C5_within_50km)
            C5_within_50km = C5_within_50km.dropna(axis=1, how='all')
            C4_sympatric = C4_sympatric._append(row, ignore_index=True)
            C5_sympatric = C5_sympatric._append(C5_within_50km, ignore_index=True)

    C4_sympatric = C4_sympatric.drop_duplicates()
    C5_sympatric = C5_sympatric.drop_duplicates()
    # print(C4_sympatric)
    C4_sympatric.to_csv(os.path.join(out_dir,'C4_sympatric2_material.txt'),sep='\t',index=False)
    C5_sympatric.to_csv(os.path.join(out_dir,'C5_sympatric2_material.txt'),sep='\t',index=False)


if __name__ == "__main__":
    # Admix_file = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Admix1_pop.list',
    #                            header=None, names=['Year', 'ID', 'Province','Lon', 'Lat', 'Alt'])
    # Admix_df = Admix_file.dropna()
    C4_file = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\C4_pop.list',
                            header=None, names=['Year', 'ID', 'Province','Lon', 'Lat', 'Alt','Ecotype'])
    C4_df = C4_file.dropna()
    C5_file = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\C5_pop.list', 
                            header=None, names=['Year', 'ID', 'Province','Lon', 'Lat', 'Alt','Ecotype'])
    C5_df = C5_file.dropna()
    out_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric'
    distance_threshold_km = 150
    main(C4_df,C5_df,distance_threshold_km,out_dir)