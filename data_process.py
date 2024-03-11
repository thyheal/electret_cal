import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from tqdm import tqdm

def check_gaussian_log(file_path):
    '''
    If valid, return True, else return False.
    '''
    try:
        with open(file_path, 'r') as file:
            log_content = file.read()
            if "Normal termination" in log_content:
                return True
            else:
                return False
    except Exception as e:
        print(f"Error: {e}")
        return False
    
def charge_calculation(dir,index):
        assert index in ['IP', 'EA'], "index should be 'IP' or 'EA'"
        neutral_path = f'{dir}/{dir}_0.log' 
        if index == 'IP':
            charge_path = f'{dir}/{dir}_p1.log'
        elif index == 'EA':
            charge_path = f'{dir}/{dir}_n1.log'
        try:
            if check_gaussian_log(charge_path) and check_gaussian_log(neutral_path):
                charge = subprocess.check_output(f"echo ` grep 'SCF Done' {charge_path} | tail -n 1 | awk '{{print $5}}' `", shell=True)
                neu = subprocess.check_output(f"echo ` grep 'SCF Done' {neutral_path} | tail -n 1 | awk '{{print $5}}' `", shell=True)
                if index == 'IP':
                    IP = (float(charge.split()[0].decode('utf-8')) - float(neu.split()[0].decode('utf-8')))* 27.2114
                    return IP
                elif index == 'EA':
                    EA = (float(neu.split()[0].decode('utf-8')) - float(charge.split()[0].decode('utf-8')))* 27.2114
                    return EA
        except:
            return None

def property_calculation(dir, index):
    assert index in ['HOMO','LUMO','HOMOn1'], "index shoule be 'HOMO', 'LUMO', 'HOMOn1'"
    if index == 'HOMOn1':
        path = f'{dir}/{dir}_n1.log'
    else:
        path = path = f'{dir}/{dir}_0.log'
    try:
        if check_gaussian_log(path):
            if index == 'LUMO':
                value = subprocess.check_output(f"echo ` grep 'virt' {path} | head -n 1 | awk '{{print $5}}' `",shell=True)

            else:
                value = subprocess.check_output(f"echo ` grep 'occ' {path} | tail -n 1 | awk '{{print $5}}' `",shell=True)
            value = (float(value.split()[0].decode('utf-8'))) * 27.2114
            return value
    except:
        return None

def prop_df(name_list,iteration,prop):
    prop_list = []
    for i in range(iteration):
        prop_list.append([])
        for j in tqdm(range(len(name_list))):
            if prop in ["IP","EA"]:
                prop_list[i].append(charge_calculation(name_list[j] + str(i),prop))
            elif prop in ["HOMO", "LUMO","HOMOn1"]:
                prop_list[i].append(property_calculation(name_list[j] + str(i),prop))
    df_prop = pd.DataFrame(prop_list, columns = name_list)
    return df_prop

def data_wash(df):
    df.replace(0, np.nan, inplace=True)
    mean_values = df.mean()
    df.fillna(mean_values, inplace=True)
    return df

def add_mean(df):
    column_means = df.mean()
    df = df.append(column_means, ignore_index=True, sort=False)
    df = df.rename(index={""})

def time_calculation(log):
    with open(log, 'r') as file:
        for line in file:
            if "Job cpu time" in line:
                cpu_time = line.split(":")[1].strip()  # 获取冒号后面的部分并去除首尾空格
                day = cpu_time.split()[0]
                hour = cpu_time.split()[2]
                min = cpu_time.split()[4]
                sec = cpu_time.split()[6]
                cpu_time = int(day) * 24 * 60 + int(hour) * 60 + int(min)  + float(sec)/60
                return cpu_time

def IP_analysis(IP_values, molecule_name):
    '''
    Return the median value of the IP values.
    '''
    ip_data = np.array(IP_values)
    # 检测和处理异常值
    ip_data_cleaned = ip_data[ip_data > 0]
    ip_data_cleaned = ip_data_cleaned[ip_data_cleaned != None]

    # 创建箱线图来展示数据分布
    plt.boxplot(ip_data_cleaned)
    plt.title(f'{molecule_name} IP distribution')
    plt.xlabel('IP')
    plt.ylim(0, 15)
    plt.ylabel('eV')
    plt.show()
    # 计算平均数和中位数
    mean_value = np.mean(ip_data_cleaned)
    median_value = np.median(ip_data_cleaned)
    min_value = np.min(ip_data_cleaned)
    max_value = np.max(ip_data_cleaned)
    
    # 创建直方图来展示数据的概率分布
    plt.hist(ip_data_cleaned, bins=10, density=True, alpha=0.6, color='b')
    plt.title(f'{molecule_name} IP distribution')
    plt.xlabel('IP')
    plt.xlim(0, 15)
    plt.ylabel('density')
    plt.show()
    print(len(ip_data_cleaned))
    print(f'This is the IP calculated for {molecule_name}')
    print(f'min:{min_value}')
    print(f'max:{max_value}')
    print(f'avg:{mean_value}')
    print(f'med:{median_value}')
    df = pd.DataFrame({'IP_Data': ip_data_cleaned})
    print(df)
    return median_value, mean_value

def corr(data1, data2):
    return np.corrcoef(np.array(data1), np.array(data2))[0,1]

def square_fig(data1, data2):
    plt.figure(figsize=(6, 6))
    plt.scatter(data1, data2)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Scatter Plot with Pearson Correlation')
    plt.legend()
    plt.show()
 