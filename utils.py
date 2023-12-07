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

def IP_calculation(dir):
    csv = f'{dir}.csv'
    open(csv, 'w').close()
    try:
        if os.path.isfile(f'{dir}/{dir}_0.log') and os.path.isfile(f'{dir}/{dir}_1.log') :
            with open(f'{dir}/{dir}_0.log') as f:
                file = f.read()
                flag1 = file.find("Normal termination")
                flag2 = file.find("Error termination")
            with open(f'{dir}/{dir}_1.log') as f:
                file = f.read()
                flag3 = file.find("Normal termination")
                flag4 = file.find("Error termination")
            if flag1 != -1 and flag3 != -1:
                os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` > {1}".format(f'{dir}/{dir}_0.log',csv))
                os.system("echo ` grep 'SCF Done' {0} | tail -n 1 | awk '{{print $5}}' ` >> {1}".format(f'{dir}/{dir}_1.log',csv))
            with open(csv, 'r') as f:
                lines = f.readlines()
                cation = lines[-2].strip()
                neutral = lines[-1].strip()
                IP = (float(cation) - float(neutral)) * 27.2114
            os.system('rm {0}'.format(csv))
            print (dir," cation energy(Ha):", cation, " neutral energy(Ha):", neutral, 'IP(eV):', IP)
            os.system("echo  {0},{1} eV >> ip_val.csv".format(dir ,IP))
    except Exception as e:
        print(f"Error: {e}")
        os.system("echo  {0},Error eV >> ip_val.csv".format(dir))

def flow_pos(molecule_name):
    xyz2gjf(xyz_name = molecule_name + "_0.xyz", gjf_header = "header_pos_nonopt.xyz")
    gjf2sh(dir = molecule_name, gjf_name = molecule_name + "_0.gjf")
    
def flow_neu(molecule_name):
    log2xyz(dir = molecule_name, log_name = molecule_name + "_0.log")
    xyz2gjf(xyz_name = molecule_name + "_1.xyz", gjf_header = "header_neu_nonopt.xyz")
    gjf2sh(dir = molecule_name, gjf_name = molecule_name + "_1.gjf")
    # log2xyz(dir = molecule_name, log_name = molecule_name + "_1.log")

def main():
    current_time = datetime.datetime.now()
    log_file_name = current_time.strftime("%Y-%m-%d_%H_%M.log")
    logging.basicConfig(filename=log_file_name, level=logging.INFO, format='%(message)s')
    logging.info('Dear user, you are using ip_cal/utils.py created by Yuhan GU.')
    smile2xyz(xyz_name = "NH3_0.xyz",smile = "N")
    flow_pos(molecule_name="NH3")

def next():
    flow_neu(molecule_name="NH3")

if __name__ == "__main__":
    for i in range(3,9):
        flow_pos(molecule_name=f"A{i}")
    for i in range(0, 4):
        flow_neu(molecule_name=f"mol{i}")
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
    