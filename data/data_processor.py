from typing import List, Dict, Optional, Union, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from tqdm import tqdm
import os
class DataProcessor:
    """数据处理核心类，负责处理高斯计算结果和数据分析。

    该类提供了一系列方法来处理和分析高斯计算的输出结果，包括：
    - 日志文件验证
    - 电荷计算（IP/EA）
    - 属性计算（HOMO/LUMO）
    - 数据清洗和统计分析
    - 结果可视化
    """

    @staticmethod
    def check_gaussian_log(file_path: Union[str, Path]) -> bool:
        """检查高斯输出日志文件是否有效。

        Args:
            file_path: 日志文件路径

        Returns:
            bool: 如果文件有效返回True，否则返回False
        """
        try:
            with open(file_path, 'r') as file:
                return "Normal termination" in file.read()
        except Exception as e:
            print(f"Error: {str(e)}")
            return False

    @staticmethod
    def charge_calculation(dir_name: str, index: str, parent_dir: Optional[str] = None) -> Optional[float]:
        """计算分子的电荷相关属性（IP或EA）。

        Args:
            dir_name: 目录名
            index: 计算类型 ('IP' 或 'EA')
            parent_dir: 父目录路径（可选）

        Returns:
            Optional[float]: 计算结果，如果计算失败返回None
        """
        if index not in ['IP', 'EA']:
            raise ValueError("index必须是'IP'或'EA'")

        # 构建正确的路径
        base_path = Path(parent_dir) / dir_name if parent_dir else Path(dir_name)
        neutral_path = base_path / f"{dir_name}_neu.log"
        charge_path = base_path / f"{dir_name}_{'pos' if index == 'IP' else 'neg'}.log"
        try:
            if not (DataProcessor.check_gaussian_log(neutral_path) and 
                    DataProcessor.check_gaussian_log(charge_path)):
                # print('abnormal')
                return None

            def get_scf_energy(path: Path) -> float:
                cmd = f"grep 'SCF Done' {path} | tail -n 1 | awk '{{print $5}}'"
                result = subprocess.check_output(cmd, shell=True)
                return float(result.decode().strip())

            neutral_energy = get_scf_energy(neutral_path)
            charge_energy = get_scf_energy(charge_path)

            # 转换为电子伏特单位
            if index == 'IP':
                return (charge_energy - neutral_energy) * 27.2114
            else:  # EA
                return (neutral_energy - charge_energy) * 27.2114

        except Exception as e:
            print(f"计算{index}时出错: {str(e)}")
            return None

    @staticmethod
    def trapdepth_calculation(dir_name: str, index: str, parent_dir: Optional[str] = None) -> Optional[float]:
        """计算分子的陷阱深度（HOMO_trapdepth或LUMO_trapdepth）。

        Args:
            dir_name: 目录名
            index: 计算类型 ('HOMO_trapdepth' 或 'LUMO_trapdepth')
            parent_dir: 父目录路径（可选）

        Returns:
            Optional[float]: 计算结果，如果计算失败返回None
        """
        if index not in ['HOMO_trapdepth', 'LUMO_trapdepth']:
            raise ValueError("index必须是'HOMO_trapdepth'或'LUMO_trapdepth'")

        try:
            if index == 'HOMO_trapdepth':
                homo = DataProcessor.property_calculation(dir_name, 'HOMO', parent_dir)
                homo_1 = DataProcessor.property_calculation(dir_name, 'HOMO-1', parent_dir)
                if homo is not None and homo_1 is not None:
                    return homo - homo_1
            else:  # LUMO_trapdepth
                lumo_1 = DataProcessor.property_calculation(dir_name, 'LUMO+1', parent_dir)
                lumo = DataProcessor.property_calculation(dir_name, 'LUMO', parent_dir)
                if lumo_1 is not None and lumo is not None:
                    return lumo_1 - lumo
            return None

        except Exception as e:
            print(f"计算{index}时出错: {str(e)}")
            return None

    @staticmethod
    def property_calculation(dir_name: str, index: str, parent_dir: Optional[str] = None) -> Optional[float]:
        """计算分子的特定属性（HOMO/HOMO-1/LUMO/LUMO+1）。

        Args:
            dir_name: 目录名
            index: 属性类型 ('HOMO', 'HOMO-1', 'LUMO', 或 'LUMO+1')
            parent_dir: 父目录路径（可选）

        Returns:
            Optional[float]: 计算结果，如果计算失败返回None
        """
        if index not in ['HOMO', 'HOMO-1', 'LUMO', 'LUMO+1']:
            raise ValueError("index必须是'HOMO'、'HOMO-1'、'LUMO'或'LUMO+1'")

        # 构建正确的路径
        base_path = Path(parent_dir) / dir_name if parent_dir else Path(dir_name)
        log_path = base_path / f"{dir_name}_neu.log"

        try:
            if not DataProcessor.check_gaussian_log(log_path):
                return None

            # 根据不同的轨道类型设置不同的提取命令
            if index in ['HOMO', 'HOMO-1']:
                # 获取最后一行occ的内容
                cmd = f"grep 'Alpha  occ. eigenvalues' {log_path} | tail -n 1"
                result = subprocess.check_output(cmd, shell=True)
                values = result.decode().strip().split()
                
                # 提取数值
                numbers = [float(val) for val in values if val.replace('.', '').replace('-', '').isdigit()]
                
                if len(numbers) == 1 and index == 'HOMO-1':
                    # 如果只有一个数值且需要HOMO-1，获取上一行的最后一个值
                    cmd = f"grep 'Alpha  occ. eigenvalues' {log_path} | tail -n 2 | head -n 1"
                    result = subprocess.check_output(cmd, shell=True)
                    prev_values = result.decode().strip().split()
                    prev_numbers = [float(val) for val in prev_values if val.replace('.', '').replace('-', '').isdigit()]
                    if prev_numbers:
                        return prev_numbers[-1] * 27.2114
                        # return prev_numbers[-1]

                    return None
                
                # HOMO是最后一个值，HOMO-1是倒数第二个值
                if len(numbers) >= (2 if index == 'HOMO-1' else 1):
                    value_index = -1 if index == 'HOMO' else -2
                    return numbers[value_index] * 27.2114
                    # return numbers[value_index]
                return None
                
            else:  # LUMO or LUMO+1
                # 获取最后一个HOMO行的下一行
                cmd = f"grep 'Alpha  occ. eigenvalues' {log_path} | tail -n 1 | awk '{{print NR}}'" # 获取最后一个HOMO行的行号
                result = subprocess.check_output(cmd, shell=True)
                homo_line_num = int(result.decode().strip())
                
                # 获取LUMO行（HOMO行的下一行）
                cmd = f"sed -n '{homo_line_num+1}p' {log_path}"
                result = subprocess.check_output(cmd, shell=True)
                values = result.decode().strip().split()
                
                # 提取数值
                numbers = [float(val) for val in values if val.replace('.', '').replace('-', '').isdigit()]
                
                if len(numbers) == 1 and index == 'LUMO+1':
                    # 如果只有一个数值且需要LUMO+1，获取下一行的值
                    cmd = f"sed -n '{homo_line_num+2}p' {log_path}"
                    result = subprocess.check_output(cmd, shell=True)
                    next_values = result.decode().strip().split()
                    next_numbers = [float(val) for val in next_values if val.replace('.', '').replace('-', '').isdigit()]
                    if next_numbers:
                        return next_numbers[0] * 27.2114
                    return None
                
                # LUMO是第一个值，LUMO+1是第二个值
                if len(numbers) >= (2 if index == 'LUMO+1' else 1):
                    value_index = 0 if index == 'LUMO' else 1
                    return numbers[value_index] * 27.2114
                return None

        except Exception as e:
            print(f"计算{index}时出错: {str(e)}")
            return None

    @staticmethod
    def create_property_dataframe(name_list: List[str], 
                                iteration: int, 
                                prop: str,
                                parent_dir: Optional[str] = None,) -> pd.DataFrame:
        """为指定的属性创建数据框。

        Args:
            name_list: 分子名称列表
            iteration: 迭代次数
            prop: 属性类型
            parent_dir: 父目录路径（可选）

        Returns:
            pd.DataFrame: 包含计算结果的数据框
        """
        prop_list = []
        for i in range(iteration):
            current_iteration = []
            for name in tqdm(name_list):
                # 构建正确的目录名，确保与文件系统中的目录结构匹配
                dir_name = f"{name}{i}"
                if prop in ["IP", "EA"]:
                    value = DataProcessor.charge_calculation(dir_name, prop, parent_dir)
                elif prop in ["HOMO_trapdepth", "LUMO_trapdepth"]:
                    value = DataProcessor.trapdepth_calculation(dir_name, prop, parent_dir)
                else:
                    value = DataProcessor.property_calculation(dir_name, prop, parent_dir)
                current_iteration.append(value)
            prop_list.append(current_iteration)

        # 创建数据框，使用分子名称作为列名，迭代序号作为行索引
        df = pd.DataFrame(prop_list, columns=name_list)
        df.index = [f"Iteration_{i}" for i in range(iteration)]
        return df

    @staticmethod
    def clean_data(df: pd.DataFrame) -> pd.DataFrame:
        """清理数据框中的异常值。

        Args:
            df: 输入数据框

        Returns:
            pd.DataFrame: 清理后的数据框
        """
        df = df.replace(0, np.nan)
        return df.fillna(df.mean())

    @staticmethod
    def analyze_property(values: List[float], 
                        molecule_name: str, 
                        property_name: str) -> Tuple[float, float]:
        """分析特定属性的统计特征并生成可视化。

        Args:
            values: 属性值列表
            molecule_name: 分子名称
            property_name: 属性名称

        Returns:
            Tuple[float, float]: (中位数, 平均值)
        """
        # 数据清理
        clean_values = np.array([v for v in values if v is not None and v > 0])

        # 创建箱线图
        plt.figure(figsize=(8, 6))
        plt.boxplot(clean_values)
        plt.title(f'{molecule_name} {property_name} Distribution')
        plt.ylabel('eV')
        plt.show()

        # 创建直方图
        plt.figure(figsize=(8, 6))
        plt.hist(clean_values, bins=10, density=True, alpha=0.6)
        plt.title(f'{molecule_name} {property_name} Distribution')
        plt.xlabel(property_name)
        plt.ylabel('Density')
        plt.show()

        # 计算统计量
        median = np.median(clean_values)
        mean = np.mean(clean_values)
        
        print(f"\n{molecule_name} {property_name} Statistics:")
        print(f"Number of valid samples: {len(clean_values)}")
        print(f"Minimum: {np.min(clean_values):.2f}")
        print(f"Maximum: {np.max(clean_values):.2f}")
        print(f"Mean: {mean:.2f}")
        print(f"Median: {median:.2f}")

        return median, mean

    @staticmethod
    def correlation_analysis(data1: List[float], 
                           data2: List[float], 
                           label1: str, 
                           label2: str) -> float:
        """计算两组数据之间的相关性并生成散点图。

        Args:
            data1: 第一组数据
            data2: 第二组数据
            label1: 第一组数据的标签
            label2: 第二组数据的标签

        Returns:
            float: 皮尔逊相关系数
        """
        correlation = np.corrcoef(data1, data2)[0, 1]

        plt.figure(figsize=(8, 8))
        plt.scatter(data1, data2, alpha=0.5)
        plt.xlabel(label1)
        plt.ylabel(label2)
        plt.title(f'Correlation Plot (r = {correlation:.3f})')
        plt.grid(True)
        plt.show()

        return correlation
    @staticmethod
    def log2xyz(log_name: str, parent_dir: Optional[str] = None, suffix: str = 'optstruc') -> bool:
        """从高斯日志文件中提取分子结构并保存为xyz格式。

        Args:
            log_name: 日志文件名
            parent_dir: 父目录路径（可选）
            suffix: xyz文件名后缀（默认为'optstruc'）

        Returns:
            bool: 如果日志文件有效返回True，否则返回False
        """
        # 获取目录名和基础名称
        dir_name = log_name.split('_')[0]
        base_name = log_name.split('.')[0].split('_')[0]
        
        # 构建完整的日志文件路径
        log_dir = os.path.join(parent_dir, dir_name) if parent_dir else dir_name
        log_path = os.path.join(log_dir, log_name)
        
        # 检查日志文件有效性
        if not DataProcessor.check_gaussian_log(log_path):
            print(f"Error: {log_path} is not a valid log file.")
            return False
            
        # 构建xyz文件名和路径
        xyz_name = f"{base_name}_{suffix}"
        temp_xyz_path = os.path.join(log_dir, f"{xyz_name}.xyz")
        
        # 使用obabel转换格式
        os.system(f'obabel {log_path} -ig09 -oxyz -O {temp_xyz_path}')
        
        # 读取并修改xyz文件
        try:
            with open(temp_xyz_path, 'r') as file:
                lines = file.readlines()
            lines[1] = '\n'  # 替换第二行为空行
            with open(temp_xyz_path, 'w') as file:
                file.writelines(lines)
        except Exception as e:
            print(f"处理xyz文件时出错: {str(e)}")
            return False
        
        # 确保目标目录存在
        final_dir = os.path.join(parent_dir, base_name) if parent_dir else base_name
        os.makedirs(final_dir, exist_ok=True)
        
        # 移动文件到最终位置
        final_path = os.path.join(final_dir, f"{xyz_name}.xyz")
        os.rename(temp_xyz_path, final_path)
        
        return True