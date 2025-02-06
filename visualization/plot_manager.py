from typing import List, Optional, Union
from pathlib import Path
from datetime import datetime
import os
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw

class PlotManager:
    """绘图管理类，负责处理分子结构和数据的可视化。

    该类提供了一系列方法来可视化分子结构和数据，包括：
    - 多分子结构的绘制
    - 属性数据的可视化
    - 图像的保存和导出

    Attributes:
        None
    """

    @staticmethod
    def generate_filename(prefix: str = "output", suffix: str = "", parent_dir: Optional[Union[str, Path]] = None, with_timestamp: bool = True) -> str:
        """生成文件名。

        Args:
            prefix: 文件名前缀
            suffix: 文件名后缀（不包含扩展名）
            parent_dir: 父目录路径（可选）
            with_timestamp: 是否包含时间戳

        Returns:
            str: 生成的文件名（包含完整路径）
        """
        filename = prefix
        if suffix:
            filename = f"{filename}_{suffix}"
        if with_timestamp:
            time_string = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            filename = f"{filename}_{time_string}"
        filename = f"{filename}.pdf"

        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
            return os.path.join(parent_dir, filename)
        return filename

    @staticmethod
    def plot_molecules(smiles_list: List[str],
                      name_list: List[str],
                      properties: Optional[List[float]] = None,
                      save_path: Optional[Union[str, Path]] = None,
                      parent_dir: Optional[Union[str, Path]] = None,
                      cols: int = 10,
                      with_timestamp: bool = True) -> None:
        """绘制多个分子结构。

        Args:
            smiles_list: SMILES字符串列表
            name_list: 分子名称列表
            properties: 分子属性值列表（可选）
            save_path: 保存路径（可选）
            parent_dir: 父目录路径（可选）
            cols: 每行显示的分子数量
            with_timestamp: 是否在文件名中包含时间戳

        Raises:
            ValueError: 如果输入参数无效
        """
        if len(smiles_list) != len(name_list):
            raise ValueError("SMILES列表和名称列表长度必须相同")
        if properties and len(properties) != len(smiles_list):
            raise ValueError("属性列表长度必须与SMILES列表相同")

        # 计算行数和列数
        num_rows = (len(smiles_list) + cols - 1) // cols
        fig, axs = plt.subplots(num_rows, cols, figsize=(cols*5, num_rows*5))

        # 确保axs是二维数组
        if num_rows == 1:
            axs = axs.reshape(1, -1)

        # 绘制分子结构
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                ax = axs[i // cols, i % cols]
                img = Draw.MolToImage(mol, size=(450, 450))
                ax.imshow(img)
                ax.axis('off')

                # 设置标题
                title = name_list[i]
                if properties:
                    title = f"{title}, {properties[i]:.2f}"
                ax.set_title(title, fontsize=30)

        # 隐藏多余的子图
        for i in range(len(smiles_list), num_rows * cols):
            row = i // cols
            col = i % cols
            axs[row, col].axis('off')

        # 调整布局
        plt.tight_layout()

        # 保存图像
        if save_path:
            if with_timestamp:
                base, ext = os.path.splitext(save_path)
                time_string = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                save_path = f"{base}_{time_string}{ext}"
            if parent_dir:
                os.makedirs(parent_dir, exist_ok=True)
                save_path = os.path.join(parent_dir, os.path.basename(save_path))
        else:
            save_path = PlotManager.generate_filename(parent_dir=parent_dir)
        plt.savefig(save_path, dpi=300, format='pdf')
        plt.close()

    @staticmethod
    def plot_property_distribution(values: List[float],
                                 property_name: str,
                                 save_path: Optional[Union[str, Path]] = None,
                                 parent_dir: Optional[Union[str, Path]] = None) -> None:
        """绘制属性分布图。

        Args:
            values: 属性值列表
            property_name: 属性名称
            save_path: 保存路径（可选）
            parent_dir: 父目录路径（可选）

        Raises:
            ValueError: 如果输入数据无效
        """
        if not values:
            raise ValueError("属性值列表不能为空")

        plt.figure(figsize=(10, 6))
        plt.hist(values, bins=30, density=True, alpha=0.7)
        plt.xlabel(property_name)
        plt.ylabel('Density')
        plt.title(f'{property_name} Distribution')
        plt.grid(True, alpha=0.3)

        # 保存图像
        save_path = save_path or PlotManager.generate_filename(f"{property_name}_dist", parent_dir=parent_dir)
        plt.savefig(save_path, dpi=300, format='pdf')
        plt.close()

    @staticmethod
    def plot_correlation(x_values: List[float],
                        y_values: List[float],
                        x_label: str,
                        y_label: str,
                        save_path: Optional[Union[str, Path]] = None,
                        parent_dir: Optional[Union[str, Path]] = None) -> None:
        """绘制相关性散点图。

        Args:
            x_values: X轴数据
            y_values: Y轴数据
            x_label: X轴标签
            y_label: Y轴标签
            save_path: 保存路径（可选）
            parent_dir: 父目录路径（可选）

        Raises:
            ValueError: 如果输入数据无效
        """
        if len(x_values) != len(y_values):
            raise ValueError("X轴和Y轴数据长度必须相同")

        plt.figure(figsize=(10, 10))
        plt.scatter(x_values, y_values, alpha=0.6)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f'{x_label} vs {y_label}')
        plt.grid(True, alpha=0.3)

        # 保存图像
        save_path = save_path or PlotManager.generate_filename(f"correlation_{x_label}_{y_label}", parent_dir=parent_dir)
        plt.savefig(save_path, dpi=300, format='pdf')
        plt.close()