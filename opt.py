# 导入必要的库
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from pathlib import Path
import time

# 添加项目根目录到系统路径
project_root = Path().absolute().parent
sys.path.append(str(project_root))

# 导入项目相关模块
from core.molecule_processor import MoleculeProcessor
from core.gaussian_calculator import GaussianCalculator
from visualization.plot_manager import PlotManager
from utils.cluster_monitor import ClusterMonitor

# some calculations for wang's paper
base_list = [
'C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)(OC(F)(F)(F))',

]
name_list = ['FFKM']
set_name = 'opt'
iteration = 3
calculator_params = {
    'parent_dir': set_name,
    'method': 'CAM-B3LYP',  # DFT泛函
    'basis': '6-31G(d,p)',  # 基组
    'opt': True,  # 进行结构优化
    'dispersion': False,  # 包含色散校正
    'polar': False,  # 不计算极化率
    'volume': False,  # 不计算体积
    'pcm': True,  # 不使用PCM溶剂化模型
    'eps': 4.9,  # PCM模型的介电常数
    'wfn': True,  # 输出波函数文件
    'debug': True  # 非调试模式
}

# 使用PlotManager绘制分子结构
PlotManager.plot_molecules(
    smiles_list=base_list,
    name_list=name_list,
    parent_dir=set_name,
    save_path='molecular_structures.pdf',
    cols=5,  # 每行显示4个分子
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
    parent_dir=set_name,
    cols=5,  # 每行显示4个分子
    with_timestamp=False
)

# ... earlier code remains the same ...

# 为每个FFKM分子生成多个xyz文件
for i, (mol_smiles, mol_name) in enumerate(zip(FFKM_list, name_list)):
    success_count = 0
    attempt = 0
    max_attempts = 20  # 最大尝试次数，防止无限循环
    
    while success_count < iteration and attempt < max_attempts:
        try:
            current_seed = 42 + attempt
            
            # 修改文件命名格式，使用 init 作为后缀
            xyz_filename = f"{mol_name}{success_count}_init.xyz"
            
            MoleculeProcessor.smile2xyz(
                xyz_name=xyz_filename,
                smile=mol_smiles,
                randomSeed=current_seed,
                parent_dir=set_name
            )
            
            success_count += 1
            # print(f"成功生成 {mol_name} 的第 {success_count} 个初始构象")
            
        except Exception as e:
            print(f"尝试生成 {mol_name} 第 {success_count + 1} 个构象时失败: {e}")
        
        attempt += 1
    
    if success_count < iteration:
        print(f"警告: {mol_name} 只成功生成了 {success_count} 个构象")

# 为每个分子准备高斯计算文件
monitor = ClusterMonitor()

for conf_id in range(iteration):
    try:
        for name in name_list:
            # while monitor.is_queue_full(4):
            #     time.sleep(60)
        # 为每个构象准备2种电荷态的计算文件
            xyz_file = f"{name}{conf_id}_init.xyz"
            # 准备阳离子态计算文件（用于IP计算）
            cation_calc = GaussianCalculator(**calculator_params, charge='pos')
            cation_calc._prepare_calculation(xyz_file)
            cation_calc._generate_input_file()
            cation_calc._generate_shell_script()
            # cation_calc._submit_job()
                    
    except Exception as e:
        print(f"处理{name}时出错: {str(e)}")
        