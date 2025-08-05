'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-08-31 14:56:03
LastEditors: Ne0tea
LastEditTime: 2025-05-01 15:25:41
'''
import pandas as pd
from geopy.distance import great_circle
import numpy as np
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point

pd.set_option('display.max_rows', None)

def calculate_geographic_center(coords):
    # 将经纬度转换为弧度
    latitudes = np.radians([float(lat) for lat, lon in coords])
    longitudes = np.radians([float(lon) for lat, lon in coords])

    x = np.cos(latitudes) * np.cos(longitudes)
    y = np.cos(latitudes) * np.sin(longitudes)
    z = np.sin(latitudes)
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    z_mean = np.mean(z)

    lon_mean = np.degrees(np.arctan2(y_mean, x_mean))
    hyp = np.sqrt(x_mean**2 + y_mean**2)
    lat_mean = np.degrees(np.arctan2(z_mean, hyp))

    return [lat_mean, lon_mean]

def within_50km(lat1, lon1, lat2, lon2,distance_threshold_km):
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    return great_circle(point1, point2).kilometers <= distance_threshold_km

def group_stats(group):
    coords = list(zip(group['Lat'], group['Lon']))
    center = calculate_geographic_center(coords)
    row_count = len(group)
    ecotype_counts = group['Ecotype'].value_counts().to_dict()
    # 创建结果字典
    result = {
        'Center_Lat': center[0],
        'Center_Lon': center[1],
        'Count': row_count
    }
    # 添加 Ecotype 计数
    for ecotype in ecotype_counts:
        result[f'{ecotype}'] = ecotype_counts[ecotype]
    return pd.Series(result)

def find_cluster(clusters,lat, lon,distance_threshold_km):
    for cluster in clusters:
        cluster_coords=[(point['Lat'],point['Lon']) for point in clusters[cluster]]
        mean_point=calculate_geographic_center(cluster_coords)
        if great_circle((lat, lon), (mean_point[0], mean_point[1])).km <= distance_threshold_km:
            return cluster
    return None

# 计算距离并划分亚群
def cluster_by_distance(df, distance_threshold_km):
    clusters = {}
    count=0
    for idx, row in df.iterrows():
        id, lat, lon = row['ID'], row['Lat'], row['Lon']
        cluster = find_cluster(clusters, lat, lon, distance_threshold_km)
        if cluster is None:
            count+=1
            clusters['Group_'+str(count)]=[row]
            continue
        clusters[cluster].append(row)
    return clusters

def assign_cluster(clusters, df, distance_threshold_km):
    out_cluster={}
    cluster_mean_point={}
    for cluster in clusters:
        cluster_coords=[(point['Lat'],point['Lon']) for point in clusters[cluster]]
        cluster_mean_point[cluster]=calculate_geographic_center(cluster_coords)

    for idx, row in df.iterrows():
        id, lat, lon = row['ID'], row['Lat'], row['Lon']
        sus_cluster=[]
        for cluster in cluster_mean_point:
            if great_circle((lat, lon), (cluster_mean_point[cluster][0], cluster_mean_point[cluster][1])).km <= distance_threshold_km:
                sus_cluster.append([cluster,great_circle((lat, lon), (cluster_mean_point[cluster][0], cluster_mean_point[cluster][1])).km])
        if sus_cluster:
            cur_cluster=sorted(sus_cluster,key=lambda x:x[1])[0][0]
            if cur_cluster not in out_cluster:
                out_cluster[cur_cluster]=[row]
            else:
                out_cluster[cur_cluster].append(row)
    return out_cluster

def main(C4_df,C5_df,distance_threshold_km,out_dir):
    All_df=C4_df._append(C5_df)
    clusters_in_C4 = cluster_by_distance(C4_df, distance_threshold_km)
    clusters_in_C5 = assign_cluster(clusters_in_C4, C5_df,distance_threshold_km)
    # print([len(clusters_in_C4[x]) for x in clusters_in_C4])
    Cluster_df=pd.DataFrame()
    for cluster in clusters_in_C4:
        for row in clusters_in_C4[cluster]:
            row['Cluster']=cluster
            Cluster_df=Cluster_df._append(row, ignore_index=True)
        for row in clusters_in_C5[cluster]:
            row['Cluster']=cluster
            Cluster_df=Cluster_df._append(row, ignore_index=True)

    out_df=pd.DataFrame()
    for cluster,df in Cluster_df.groupby('Cluster'):
        cur_row=group_stats(df)
        cur_row['Cluster']=cluster
        out_df=out_df._append(cur_row, ignore_index=True)
    out_df.fillna(0,inplace=True)
    out_df=out_df.set_index('Cluster')
    out_df=out_df[(out_df['Count'] >5) & (out_df['C4']>=2)]

    result_group=out_df.index.to_list()
    Cluster_df=Cluster_df[Cluster_df['Cluster'].isin(result_group)]
    print(Cluster_df)

    material_overlap = []
    for _, ma_row in Cluster_df.iterrows():
        cur_lon,  cur_lat = ma_row['Lon'], ma_row['Lat']
        if ma_row['Ecotype'] == 'C4': continue
        for idx, cluster in out_df.iterrows():
            if idx == ma_row['Cluster']: continue
            cluster_lat, cluster_lon = cluster['Center_Lat'], cluster['Center_Lon']
            if within_50km(cur_lat, cur_lon, cluster_lat, cluster_lon, distance_threshold_km):
                material_overlap.append(ma_row['ID'])
                # print(ma_row['ID'])

    Cluster_df = Cluster_df[~Cluster_df['ID'].isin(material_overlap)]
    # print(len(Cluster_df))
    # print(Cluster_df.groupby(['Cluster','Ecotype']).size())
    out_df.to_csv(os.path.join(out_dir,'Sympatric_group_center2.txt'),sep='\t')
    Cluster_df.to_csv(os.path.join(out_dir,'Sympatric_material2.txt'),sep='\t',index=False)

    # centroids = Cluster_df.groupby('Cluster').apply(group_stats).reset_index()
    # print(centroids)
    # print(Cluster_df['Cluster'].value_counts())
    if 0:
        out_df=out_df.reset_index()
        geometry = [Point(xy) for xy in zip(out_df['Center_Lon'], out_df['Center_Lat'])]
        # 创建 GeoDataFrame
        gdf = gpd.GeoDataFrame(out_df, geometry=geometry)
        # 获取世界地图数据
        world = gpd.read_file(r"E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\ne_110m_admin_0_countries.shp")
        # print(world)
        # 只保留中国部分
        china = world[world['SOVEREIGNT'] == 'China']
        # 绘制地图
        fig, ax = plt.subplots(figsize=(10, 4))
        # 绘制中国地图
        china.plot(ax=ax, color='white', edgecolor='black')
        # 在地图上绘制点
        gdf.plot(ax=ax, color='black', markersize=100)
        # 限制显示范围，只显示中国部分
        plt.xlim(100.0, 135.0)
        plt.ylim(32.0, 45.0)
        # 添加标签
        for x, y, label in zip(gdf.geometry.x, gdf.geometry.y, gdf['Cluster']):
            ax.text(x, y, label, fontsize=12, ha='right')

        # 设置标题和轴标签
        plt.title('Geographical Points in China')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        # 显示地图
        plt.show()


if __name__ == "__main__":
    C4_sym_df = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\C4_sympatric2_material.txt', header=0)
    C5_sym_df = pd.read_table(r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric\C5_sympatric2_material.txt', header=0)
    out_dir = r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\03pop_structure\Introgression\Sympatric'
    distance_threshold_km = 120
    main(C4_sym_df,C5_sym_df,distance_threshold_km,out_dir)