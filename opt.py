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
    'wfn': True,  # 输出波函数文件
    'debug': True  # 非调试模式
}

    'opt': True,  # 进行结构优化
PlotManager.plot_molecules(
    smiles_list=base_list,
    name_list=name_list,
    'pcm': True,  # 不使用PCM溶剂化模型
    save_path='molecular_structures.pdf',
    cols=5,  # 每行显示4个分子
    'debug': True  # 非调试模式
)
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
        