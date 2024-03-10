import os
from rdkit import Chem
from rdkit.Chem import AllChem
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
def smile2xyz(xyz_name,smile,randomSeed = None):
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
    if randomSeed is not None:
        AllChem.EmbedMolecule(mol, randomSeed=randomSeed)
    else:
        AllChem.EmbedMolecule(mol)
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

def i8cpu_running():
    output = subprocess.check_output("squeue", shell=True, text=True)
    count = output.count('i8cpu')
    if count > 0:
        return True
    else:
        return False
  
