'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-08-30 14:25:15
LastEditors: Ne0tea
LastEditTime: 2024-09-05 19:46:02
'''
import sys
import os
import multiprocessing
from datetime import datetime
backup_stairpainter=r'/public/home/huangyj/software/stairway_plot_v2_backup/stairway_plot_v2.1.1/stairway_plot_es/Stairpainter.class'
# 定义读取 shell 脚本文件并提取命令的函数
def read_shell_script(file_path):
    commands_create_theta = [] #create .addTheta files
    commands_move = [] #comamnd "mv"
    commands_paint = [] #visulize .addTheta files
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#"):  # 忽略注释和空行
                if 'Stairway_unfold_training_testing7' in line or 'Stairway_fold_training_testing7' in line:
                    commands_create_theta.append(line)
                elif 'mv -f' in line:
                    commands_move.append(line)
                elif 'Stairpainter' in line or 'plot.sh' in line:
                    commands_paint.append(line)
    return commands_create_theta, commands_move, commands_paint

# 定义执行单个命令的函数
def run_command(command):
    print(f"Now command {command}")
    try:
        os.system(command)
    except:
        print(f"Error occurred while executing command: {command}")
        os.sysexit(1)

# 使用 multiprocessing.Pool 进行多进程执行
def main(sh_file,multi_thread):
    # 输入 shell 脚本文件路径
    Task_start = datetime.now()
    shell_script_path = sh_file

    commands_step1,commands_step2,commands_step3 = read_shell_script(shell_script_path)
    num_processes = int(multi_thread)

    # 使用多进程执行命令
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(run_command, commands_step1)
        pool.close()
        pool.join()

    step1_time=datetime.now()-Task_start
    last_time=datetime.now()
    print(f"step1 time:{step1_time}")

    for i in commands_step2:
        os.system(i)
    step2_time=datetime.now()-last_time
    last_time=datetime.now()
    print(f"step2 time:{step2_time}")

    for i in commands_step3:
        if 'plot.sh' in i and 'bash' in i:
            plot_file=i.strip().split(' ')[-1]
            sed_command='sed -i "s/Xmx4g/Xmx50g/g" '+plot_file
            os.system(sed_command)
        elif 'Stairpainter' in i:
            os.system(i)
            cp_command='cp '+backup_stairpainter+' '+i.strip().split(' ')[-3]
            os.system(cp_command)
            continue
        os.system(i)
    step3_time=datetime.now()-last_time
    print(f"step3 time:{step3_time}")

    print('Job finish!')

if __name__ == "__main__":
    sh_file=sys.argv[1]
    multi_thread=sys.argv[2]
    main(sh_file,multi_thread)


