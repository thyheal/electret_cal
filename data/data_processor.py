from typing import List, Dict, Optional, Union, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from tqdm import tqdm

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
                print('abnormal')
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
    def property_calculation(dir_name: str, index: str, parent_dir: Optional[str] = None) -> Optional[float]:
        """计算分子的特定属性（HOMO/LUMO/HOMOn1）。

        Args:
            dir_name: 目录名
            index: 属性类型 ('HOMO', 'LUMO', 或 'HOMOn1')
            parent_dir: 父目录路径（可选）

        Returns:
            Optional[float]: 计算结果，如果计算失败返回None
        """
        if index not in ['HOMO', 'LUMO', 'HOMOn1']:
            raise ValueError("index必须是'HOMO'、'LUMO'或'HOMOn1'")

        # 构建正确的路径
        base_path = Path(parent_dir) / dir_name if parent_dir else Path(dir_name)
        log_path = base_path / f"{dir_name}_{'neg' if index == 'HOMOn1' else 'neu'}.log"

        try:
            if not DataProcessor.check_gaussian_log(log_path):
                return None

            cmd = f"grep '{'virt' if index == 'LUMO' else 'occ'}' {log_path} | {'head' if index == 'LUMO' else 'tail'} -n 1 | awk '{{print $5}}'"
            result = subprocess.check_output(cmd, shell=True)
            return float(result.decode().strip()) * 27.2114

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

        Returns:
            pd.DataFrame: 包含计算结果的数据框
        """
        prop_list = []
        for i in range(iteration):
            current_iteration = []
            for name in tqdm(name_list):
                name_with_iteration = name + str(i)
                if prop in ["IP", "EA"]:
                    value = DataProcessor.charge_calculation(name_with_iteration, prop, parent_dir)
                else:
                    value = DataProcessor.property_calculation(name_with_iteration, prop, parent_dir)
                current_iteration.append(value)
            prop_list.append(current_iteration)

        return pd.DataFrame(prop_list, columns=name_list)

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