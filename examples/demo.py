# 导入必要的库
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

# 添加项目根目录到系统路径
project_root = Path().absolute().parent
sys.path.append(str(project_root))

# 导入项目相关模块
from core.molecule_processor import MoleculeProcessor
from core.gaussian_calculator import GaussianCalculator
from visualization.plot_manager import PlotManager

# some calculations for wang's paper
base_list = [
'C1(=O)N(CC=C)C2C(N(CC=C)C(=O)N2(CC=C))N1(CC=C)',
'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CC=C)',

'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(C)',
'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CC)',
'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CCC)',
'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CCCC)',
'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CCCCC)',

'C1(=O)N(CC=C)C(=O)N(CC=C)C(=O)N1(CCN2C(=O)N(CC=C)C(=O)N(CC=C)C2(=O))',

]
name_list = ['finalTAG','finalTAIC','finalLDAIC1','finalLDAIC2','finalLDAIC3','finalLDAIC4','finalLDAIC5','finalDD2']

# 绘制原始分子结构
PlotManager.plot_molecules(
    smiles_list=base_list,
    name_list=name_list,
    save_path='original_molecules.pdf',
    cols=4,  # 每行显示4个分子
    with_timestamp=False
)

# FFKM分子构建和交联反应参数
direction = 'L'
m = 1
n = 0
c = 3

# 进行FFKM分子构建和交联反应
FFKM_list = []
for base_mol in base_list:
    # 构建FFKM分子
    ffkm, connection_point = MoleculeProcessor.build_ffkm(
        m=m,
        n=n,
        c=c,
        connection_direction=direction
    )
    
    # 模拟交联反应
    crossed_mol = MoleculeProcessor.simulate_crosslink_reaction(
        base=base_mol,
        ffkm=ffkm,
        connection_point=connection_point
    )
    
    FFKM_list.append(crossed_mol)

# 绘制交联反应后的分子结构
PlotManager.plot_molecules(
    smiles_list=FFKM_list,
    name_list=name_list,
    save_path='crosslinked_molecules.pdf',
    cols=4,  # 每行显示4个分子
    with_timestamp=False
)

# 为每个FFKM分子生成多个xyz文件
for i, (mol_smiles, mol_name) in enumerate(zip(FFKM_list, name_list)):
    success_count = 0
    attempt = 0
    max_attempts = 20  # 最大尝试次数，防止无限循环
    
    while success_count < 3 and attempt < max_attempts:
        try:
            # 使用尝试次数作为随机种子的一部分
            current_seed = 42 + attempt
            
            # 生成xyz文件名，使用分子名称和成功次数编号
            xyz_filename = f"{mol_name}_{success_count}.xyz"
            
            # 使用MoleculeProcessor的smile2xyz方法生成xyz文件
            MoleculeProcessor.smile2xyz(
                xyz_name=xyz_filename,
                smile=mol_smiles,
                randomSeed=current_seed
            )
            success_count += 1
        except Exception as e:
            print(f"生成{mol_name}的第{success_count}个构象时出错: {str(e)}")
        attempt += 1

# 设置高斯计算参数
calculator_params = {
    'method': 'B3LYP',  # DFT泛函
    'basis': '6-31G(d)',  # 基组
    'opt': True,  # 进行结构优化
    'dispersion': True,  # 包含色散校正
    'polar': False,  # 不计算极化率
    'volume': False,  # 不计算体积
    'pcm': False,  # 不使用PCM溶剂化模型
    'eps': 0,  # PCM模型的介电常数
    'wfn': True,  # 输出波函数文件
    'debug': False  # 非调试模式
}

# 为每个分子准备高斯计算文件
for name in name_list:
    try:
        # 创建分子目录
        mol_dir = Path('examples') / name
        mol_dir.mkdir(parents=True, exist_ok=True)
        
        # 为每个构象准备三种电荷态的计算文件
        for conf_id in range(3):
            xyz_file = f"{name}_{conf_id}.xyz"
            
            # 准备中性态计算文件
            neutral_calc = GaussianCalculator(**calculator_params, charge='neu')
            neutral_calc._prepare_calculation(xyz_file)
            neutral_calc._generate_input_file()
            neutral_calc._generate_shell_script()
            
            # 准备阳离子态计算文件（用于IP计算）
            cation_calc = GaussianCalculator(**calculator_params, charge='pos')
            cation_calc._prepare_calculation(xyz_file)
            cation_calc._generate_input_file()
            cation_calc._generate_shell_script()
            
        print(f"成功为{name}准备所有计算文件")
        
    except Exception as e:
        print(f"处理{name}时出错: {str(e)}")