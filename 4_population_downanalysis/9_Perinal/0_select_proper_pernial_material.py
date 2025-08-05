'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-01-29 20:36:04
LastEditors: Ne0tea
LastEditTime: 2025-04-12 15:48:55
'''
from geopy.distance import geodesic

def filter_points(series1, series2, threshold=10):
    filtered_series1 = []

    for i in series1:
        point1=(list(i.values())[0][1],list(i.values())[0][0])
        # 判断在series2中是否存在距离在threshold以内的点
        if any(geodesic(point1, (list(x.values())[0][1],list(x.values())[0][0])).km < threshold for x in series2):
            filtered_series1.append(i)

    filtered_series2 = []
    for i in series2:
        point2=(list(i.values())[0][1],list(i.values())[0][0])
        if any(geodesic(point2, (list(x.values())[0][1],list(x.values())[0][0])).km < threshold for x in series1):
            filtered_series2.append(i)

    return filtered_series1, filtered_series2

def main(file1,file2,outfile,cutoff):
    loc1_list=[]
    with open(file1,'r') as f1:
        for i in f1:
            line=i.strip().split()
            loc1_list.append({str(line[0]):(float(line[1]),float(line[2]))})
    loc2_list=[]
    with open(file2,'r') as f2:
        for i in f2:
            line=i.strip().split()
            loc2_list.append({str(line[0]):(float(line[1]),float(line[2]))})

    filtered_series1, filtered_series2 = filter_points(loc1_list, loc2_list,cutoff)

    print("保留在series1中的点:", len(loc1_list),len(filtered_series1))
    print("保留在series2中的点:", len(loc2_list),len(filtered_series2))
    with open(outfile,'w') as of:
        of.write('id\tlong\tlat\tyear'+'\n')
        for i in filtered_series1:
            id=list(i.keys())
            loc=list(list(i.values())[0])
            id.extend(loc)
            xx=[str(x) for x in id]
            xx.append('2013')
            line='\t'.join(xx)
            of.write(line+'\n')
        for i in filtered_series2:
            id=list(i.keys())
            loc=list(list(i.values())[0])
            id.extend(loc)
            xx=[str(x) for x in id]
            xx.append('2023')
            line='\t'.join(xx)
            of.write(line+'\n')

if __name__ == "__main__":
    # loc_file=sys.argv[1]
    loc1_file=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Material_2013.txt'
    loc2_file=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Material_2023.txt'
    outfile=r'E:\Bio_analysis\Weed_genome_project\Digitaria\Population\material\Population_perinal.txt'
    cutoff=100
    main(loc1_file,loc2_file,outfile,cutoff)