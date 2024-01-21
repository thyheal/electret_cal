import shutil
from openbabel import pybel
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from datetime import datetime
import csv
import logging
import datetime
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import subprocess

'''
This is some tools used to calculate the ip of molecules.
    1.The flow of the calculation can be illustrated as follows:
      smiles -> xyz0  -> gjf0 -> sh0 -> log0 -> xyz1 -> gjf1 -> sh1 -> log1 -> ip
      some modules are used to generate the files needed in the calculation.
      such as smiles2xyz, xyz2gjf, gjf2sh, log2xyz, log2energy
    2.The calculation process should be seprarated into two parts:
        1) the first part is to calculate the ip in positive charge state.
        2) the second part is to calculate the ip in neutral charge state.
        Hence the flow was done with function flow_pos and flow_neu.

'''
# ###########################################################
def smile2xyz(xyz_name,smile):
    '''
    xyz_name: xyz file including postfixed .xyz
    smiles: smiles string
    for instance smile2xyz('CH4_0.xyz','C')
    create a CH4 folder and => dir_name
    create a CH4_0.xyz file named CH4.xyz stored in  CH4 folder
    !!! xyz file should be ended with _n.xyz, n is a index, should be 1/2.
    '''
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed = 1)
    AllChem.MMFFOptimizeMolecule(mol)
    # Draw.MolToImage(mol,size=(100,100))
    xyz_str = Chem.MolToXYZBlock(mol)
    dir_name = xyz_name.split('_')[0]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    file_path = os.path.join(dir_name, xyz_name)
    with open(file_path, "w") as file:
        file.write(xyz_str)

def xyzcheck(xyz_name, Canonsmile):
    '''
    This is used to check whether the xyz file is valid.
    Valid means the xyz file can be converted to the same smile string as input in GeoMol.
    '''
    molecule_name =  xyz_name.split('_')[0]
    os.system('obabel {2}/{0} -ixyz --osmi -O {2}/{1}.smi'.format(xyz_name,molecule_name,molecule_name))
    with open(f'{molecule_name}/{molecule_name}.smi', 'r') as file:
        lines = file.readlines()
    smile = lines[0].split('\t')[0]
    try:
        smile1 = smile.replace('@','')
        # smile1 = smile.replace('[','')
        # smile1 = smile.replace(']','')
        smile1 = Chem.CanonSmiles(smile1)
        Canonsmile1 = Canonsmile.replace('@','')
        Canonsmile1 = Chem.CanonSmiles(Canonsmile1)
        if smile1 != Canonsmile1:
            return False
        else:
            return True
    except Exception as e:
        print(f"Error in xyzcheck: {e}")
        return False

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

def log2xyz(log_name):
    '''
    If valid, return True, else return False.
    sometimes the Gaussian calculation will fail, so we need to check the log file.
    '''
    dir = log_name.split('_')[0]
    if not check_gaussian_log(os.path.join(dir, log_name)):
        print(f"Error: {dir}/{log_name} is not a valid log file.")        
    id = int(log_name.split('.')[0].split('_')[1]) + 1
    xyz_name = log_name.split('.')[0].split('_')[0] + '_' + str(id)
    os.system('obabel {0}/{1} -ig09 -oxyz -O {0}/{2}.xyz'.format(dir,log_name,xyz_name))
    with open(f'{dir}/{xyz_name}.xyz', 'r') as file:
        lines = file.readlines()
    # replace the second line with blank line
    print(lines[1])
    lines[1] = '\n'
    # rewrite the xyz file
    with open(f'{dir}/{xyz_name}.xyz', 'w') as file:
        file.writelines(lines)
    return check_gaussian_log(os.path.join(dir, log_name))

def IP_calculation(dir):
    csv = f'{dir}.csv'
    open(csv, 'w').close()
    try:
        path1 = f'{dir}/{dir}_p1.log'
        path2 = f'{dir}/{dir}_0.log'
        if check_gaussian_log(path1) and check_gaussian_log(path2):
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` > {1}".format(path1,csv))
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` >> {1}".format(path2,csv))
            with open(csv, 'r') as f:
                lines = f.readlines()
                cation = lines[-2].strip()
                neutral = lines[-1].strip()
                IP = (float(cation) - float(neutral)) * 27.2114
                os.system('rm {0}'.format(csv))
                os.system("echo  {0},{1} eV >> ip_val.csv".format(dir ,IP))
            # print (dir," cation energy(Ha):", cation, " neutral energy(Ha):", neutral, 'IP(eV):', IP)
                return(IP)
        else:
            os.system('rm {0}'.format(csv))
            return(0)
    except Exception as e:
        print(f"Error: {e}")
        os.system("echo  {0},Error eV >> ip_val.csv".format(dir))
        os.system('rm {0}'.format(csv))
        return(0)
def EA_calculation(dir):
    csv = f'{dir}.csv'
    open(csv, 'w').close()
    try:
        path1 = f'{dir}/{dir}_0.log'
        path2 = f'{dir}/{dir}_n1.log'
        if check_gaussian_log(path1) and check_gaussian_log(path2):
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` > {1}".format(path1,csv))
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` >> {1}".format(path2,csv))
            with open(csv, 'r') as f:
                lines = f.readlines()
                cation = lines[-2].strip()
                neutral = lines[-1].strip()
                IP = (float(cation) - float(neutral)) * 27.2114
                os.system('rm {0}'.format(csv))
                os.system("echo  {0},{1} eV >> ea_val.csv".format(dir ,IP))
            # print (dir," cation energy(Ha):", cation, " neutral energy(Ha):", neutral, 'IP(eV):', IP)
                return(IP)
        else:
            os.system('rm {0}'.format(csv))
            return(0)
    except Exception as e:
        print(f"Error: {e}")
        os.system("echo  {0},Error eV >> ea_val.csv".format(dir))
        os.system('rm {0}'.format(csv))
        return(0)

def time_calculation(log):
    with open(log, 'r') as file:
        for line in file:
            if "Job cpu time" in line:
                cpu_time = line.split(":")[1].strip()  # 获取冒号后面的部分并去除首尾空格
                print(f"Job cpu time: {cpu_time}")
                break

def i8cpu_running():
    output = subprocess.check_output("squeue", shell=True, text=True)
    count = output.count('i8cpu')
    if count > 0:
        return True
    else:
        return False

    
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
 