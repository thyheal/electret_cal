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
    Draw.MolToImage(mol,size=(100,100))
    xyz_str = Chem.MolToXYZBlock(mol)
    dir_name = xyz_name.split('_')[0]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    file_path = os.path.join(dir_name, xyz_name)
    with open(file_path, "w") as file:
        file.write(xyz_str)

def xyzcheck(molecule_name, Canonsmile):
    '''
    This is used to check whether the xyz file is valid.
    Valid means the xyz file can be converted to the same smile string as input in GeoMol.
    '''
    xyz_name = molecule_name + "_0.xyz"
    os.system('obabel {2}/{0} -ixyz --osmi -O {2}/{1}.smi'.format(xyz_name,molecule_name,molecule_name))
    with open(f'{molecule_name}/{molecule_name}.smi', 'r') as file:
        lines = file.readlines()
    smile = lines[0].split('\t')[0]
    try:
        smile1 = smile.replace('@','')
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
    
def xyz2gjf(xyz_name, gjf_header):
    '''
    xyz_name: xyz file including postfixed .xyz
    gjf_header: gjf header file including postfixed .xyz
    '''
    dir_name = xyz_name.split('_')[0]
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    xyz_path = os.path.join(dir_name, xyz_name)
    header_path = os.path.join('header', gjf_header)
    gjf_path = os.path.join(dir_name, xyz_name.replace(".xyz", ".gjf"))

    with open(header_path, 'r') as header_file:
        header_lines = header_file.readlines()
    with open(xyz_path, 'r') as xyz_file:
        atomic_coordinates = xyz_file.readlines()[1:]  # skip the first line
    wfname = xyz_name.replace(".xyz", ".wfn")

    with open(gjf_path, 'w') as gjf_file:
        gjf_file.writelines(header_lines)
        gjf_file.writelines(atomic_coordinates)
        gjf_file.write("\n")  # add a blank line
        gjf_file.write(f"{dir_name}/{wfname}\n") #add a line of wfname
        gjf_file.write("\n" * 4)  # add 4 blank lines, this is important for gaussian calculation

    #change the molecule name in gjf file __change__
    molename = xyz_name.split('.')[0]
    with open(gjf_path, 'r') as gjf_file:
        file_content = gjf_file.read()
        modified_content = file_content.replace("__change__", molename)
        modified_content = modified_content.replace("__dir__", dir_name)
    #write the modified content back to the file
    with open(gjf_path, 'w') as gjf_file:
        gjf_file.write(modified_content)
    state = gjf_header.split('_')[1].split('.')[0]
    # log.info(f"gjf file with the state of {state} for {dir_name}/{xyz_name} has been generated.")
    # print(f"gjf file with the state of {state} for {dir_name}/{xyz_name} has been generated.")

def gjf2sh(dir, gjf_name):
    sh_path = os.path.join(dir, gjf_name.replace(".gjf", ".sh"))
    sh_header_path = os.path.join('header', "sub_g09.xyz")
    with open(sh_header_path, 'r') as file:
        file_content = file.read()
        modified_content = file_content.replace("__dir__", dir)
        modified_content = modified_content.replace("__name__", gjf_name.replace(".gjf", ""))
        modified_content = modified_content.replace("__changegjf__", gjf_name)
        modified_content = modified_content.replace("__changelog__", gjf_name.replace(".gjf", ".log"))
    with open(sh_path, 'w') as file:
        file.write(modified_content)
    os.system("pjsub {}".format(sh_path))


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

def log2xyz(dir, log_name):
    '''
    If valid, return True, else return False.
    sometimes the Gaussian calculation will fail, so we need to check the log file.
    '''
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


def flow_pos(molecule_name):
    xyz2gjf(xyz_name = molecule_name + "_0.xyz", gjf_header = "header_pos_nonopt.xyz")
    gjf2sh(dir = molecule_name, gjf_name = molecule_name + "_0.gjf")
    
def flow_neg(molecule_name):
    xyz2gjf(xyz_name = molecule_name + "_0.xyz", gjf_header = "header_neg_nonopt.xyz")
    gjf2sh(dir = molecule_name, gjf_name = molecule_name + "_0.gjf")

def flow_neu(molecule_name):
    log2xyz(dir = molecule_name, log_name = molecule_name + "_0.log")
    xyz2gjf(xyz_name = molecule_name + "_1.xyz", gjf_header = "header_neu_nonopt.xyz")
    gjf2sh(dir = molecule_name, gjf_name = molecule_name + "_1.gjf")
    # log2xyz(dir = molecule_name, log_name = molecule_name + "_1.log")


def IP_calculation(dir):
    csv = f'{dir}.csv'
    open(csv, 'w').close()
    try:
        path1 = f'{dir}/{dir}_0.log'
        path2 = f'{dir}/{dir}_1.log'
        if check_gaussian_log(path1) and check_gaussian_log(path2):
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` > {1}".format(f'{dir}/{dir}_0.log',csv))
            os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` >> {1}".format(f'{dir}/{dir}_1.log',csv))
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


# def main():
#     current_time = datetime.datetime.now()
#     log_file_name = current_time.strftime("%Y-%m-%d_%H_%M.log")
#     logging.basicConfig(filename=log_file_name, level=logging.INFO, format='%(message)s')
#     logging.info('Dear user, you are using ip_cal/utils.py created by Yuhan GU.')
#     smile2xyz(xyz_name = "NH3_0.xyz",smile = "N")
#     flow_pos(molecule_name="NH3")


# if __name__ == "__main__":
#     for i in range(3,9):
#         flow_pos(molecule_name=f"A{i}")
#     for i in range(0, 4):
#         flow_neu(molecule_name=f"mol{i}")
    # for i in range(0, 10):
    #     IP_calculation(dir=f"mol{i}")
# if __name__ == "__main__":
#     smile2xyz(xyz_name = "NH3_0.xyz",smile = "N")
#     xyz2gjf(xyz_name = "NH3_0.xyz", gjf_header = "header_pos_nonopt.xyz")
#     gjf2sh(dir = "NH3", gjf_name = "NH3_0.gjf")
#     log2xyz(dir = "NH3", log_name = "NH3_0.log")
#     xyz2gjf(xyz_name = "NH3_1.xyz", gjf_header = "header_neu_nonopt.xyz")
#     gjf2sh(dir = "NH3", gjf_name = "NH3_1.gjf")
#     log2xyz(dir = "NH3", log_name = "NH3_1.log")
#     IP_calculation(dir = "NH3")
    
# name_list = ['EHOPA','AEPY','AEP','APN','DBE','DIPEDA','OA','mXD','S','A','TAEA']
name_list = ['EHOPA','AEPY','AEP','APN','DBE','DIPEDA','OA','mXD','S','A']
smiles_list = [
'CCCCC(CC)COCCCNC(C(F)(F)C1(F)OC(F)(F)C(F)(F)C1(C(F)(C(O)=O)F)F)=O',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCC2=CN=CC=C2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN2CCCCC2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NC2=CC(C#N)=CC(C#N)=C2)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN(CCCC)CCCC)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCN(C(C)C)C(C)C)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCCCCCCCC)=O)F)O1',
'FC1(F)C(F)(F)C(F)(C(F)(C(O)=O)F)C(F)(C(F)(C(NCC2=CC(CN)=CC=C2)=O)F)O1',
'FC1(C(F)(C(F)(F)F)F)C(OC(F)(F)C(F)1F)(C(F)(C(F)(F)F)F)F',
'FC1(C(F)(C(F)(F)F)F)C(OC(F)(F)C(F)1F)(C(F)(C(O)=O)F)F',
]

smiles_dict = dict(zip(name_list, smiles_list))


names = locals()
IP_PCM = {  'EHOPA':7.38,
            'AEPY':7.63,
            'AEP':6.18,
            'APN':8.18,
            'DBE':6.27,
            'DIPEDA':6.21,
            'OA':8.15,
            'mXD':6.78,
            'S':11.32,
            'A':9.73,
        }
P_sp = {    'EHOPA':0.686,
            'AEPY':0.855,
            'AEP':0.862,
            'APN':0.558,
            'DBE':0.746,
            'DIPEDA':0.715,
            'OA':0.576,
            'mXD':0.959,
            'S':0.135,
            'A':0.211,
        }
N_sp = {    'EHOPA':0.5,
            'AEPY':0.145,
            'AEP':0.815,
            'APN':0.409,
            'DBE':0.676,
            'DIPEDA':0.592,
            'OA':0.481,
            'mXD':0.954,
            'S':0.088,
            'A':0.237,
        }

#IP calculation gaussian
for _ in name_list:
    count = 0
    for i in range(20,40):
        if xyzcheck(f'{_}{i}', smiles_dict[_]):
            flow_pos(molecule_name=f'{_}{i}')
            time.sleep(150)
            if log2xyz(dir=f'{_}{i}', log_name=f'{_}{i}_0.log'):
                flow_neu(molecule_name=f'{_}{i}')
                time.sleep(150)
                count+=1
            else:
                continue
            # print(f'{_}{i} is a valid smile string.')
        else:
            continue
            # print(f'{_}{i} is not a valid smile string.')
    print(_, count)

# IP calculation
def process_data(prefix, num):
    IP_list = []
    count_list = []
    try:
        for i in range(1, num):
            identifier = f'{prefix}{i}'
            if xyzcheck(identifier,smiles_dict[prefix]):
                IP_list.append(IP_calculation(identifier))
                count_list.append(i)
                print(f'{identifier} is ok')
            else:
                print(f'{identifier} is error')
    except Exception as e:
        print(f"Error: {e}")
    return IP_list

def corr(data1, data2):
    return np.corrcoef(np.array(data1), np.array(data2))[0,1]

names = locals()
for _ in name_list:
    names[f'IP_{_}'] = process_data(_, 20)

# IP analysis graph
IP_mean,IP_med = {},{}
for _ in name_list:
    names[f'IP_med_{_}'], names[f'IP_avg_{_}']= IP_analysis(names[f'IP_{_}'],_)
    IP_mean[_] = names[f'IP_avg_{_}']
    IP_med[_] = names[f'IP_med_{_}']

def square_fig(data1, data2):
    plt.figure(figsize=(6, 6))
    plt.scatter(data1, data2)
    # 添加标签和标题
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Scatter Plot with Pearson Correlation')
    # 添加图例
    plt.legend()
    # 显示图表
    plt.show()

psp = []
for _ in name_list:
    psp.append(P_sp[_])
ippcm = []
for _ in name_list:
    ippcm.append(IP_PCM[_])
ipmean = []
for _ in name_list:
    ipmean.append(IP_mean[_])
ipmed = []
for _ in name_list:
    ipmed.append(IP_med[_])
nsp = []
for _ in name_list:
    nsp.append(N_sp[_])


for _ in name_list:
    names[f'valid_{_}'] = 0
for _ in name_list:
    for i in range(0,40):
        if xyzcheck(f'{_}{i}', smiles_dict[_]):
            names[f'valid_{_}'] += 1

for _ in name_list:
    print(f'{_} valid: {names[f"valid_{_}"]}')


name_list = ['DIPEDA','OA','mXD','S','A']
for _ in name_list:
    count = 0
    for i in range(20,40):
        if xyzcheck(f'{_}{i}', smiles_dict[_]):
            flow_pos(molecule_name=f'{_}{i}')
            time.sleep(150)
            if log2xyz(dir=f'{_}{i}', log_name=f'{_}{i}_0.log'):
                flow_neu(molecule_name=f'{_}{i}')
                time.sleep(100)
                count+=1
            else:
                continue
            # print(f'{_}{i} is a valid smile string.')
        else:
            continue
            # print(f'{_}{i} is not a valid smile string.')
    print(_, count)