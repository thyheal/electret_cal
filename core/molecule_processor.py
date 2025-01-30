from typing import Tuple, List
from rdkit import Chem
from rdkit.Chem import AllChem
import os

class MoleculeProcessor:
    """分子结构处理核心类，负责生成和处理分子结构。

    该类提供了一系列方法来处理分子结构，包括：
    - FFKM分子的构建
    - 交联反应的模拟
    - 分子结构的验证和优化

    Attributes:
        None
    """

    @staticmethod
    def build_ffkm(m: int, n: int, c: int, connection_direction: str) -> Tuple[str, int]:
        """构建FFKM分子结构。

        Args:
            m: 含氧基团数量
            n: 不含氧基团数量
            c: 连接碳原子数量
            connection_direction: 连接方向 ('L' 或 'R')

        Returns:
            Tuple[str, int]: (FFKM分子SMILES字符串, 连接点位置)

        Raises:
            ValueError: 如果连接方向无效
        """
        if connection_direction not in ['L', 'R']:
            raise ValueError("连接方向必须是'L'或'R'")

        FFKM = n * 'C(F)(F)' + m * 'C(F)(F)C(OC(F)(F)(F))(F)'

        if connection_direction == 'L':
            connection_point = 0
            FFKM = c * 'C' + FFKM + 'F'
        else:  # R
            FFKM = 'F' + FFKM + 'C' * c
            if c == 0:
                connection_point = 3 * n + 10 * (m - 1) + 3 + 1
            else:
                connection_point = 3 * n + 10 * m + c

        return FFKM, connection_point

    @staticmethod
    def simulate_crosslink_reaction(base: str, ffkm: str, connection_point: int) -> str:
        """模拟交联反应过程。

        Args:
            base: 基础分子的SMILES字符串
            ffkm: FFKM分子的SMILES字符串
            connection_point: 连接点位置

        Returns:
            str: 反应产物的SMILES字符串

        Raises:
            ValueError: 如果分子结构无效
        """
        try:
            mol = Chem.MolFromSmiles(base)
            if mol is None:
                raise ValueError("无效的基础分子结构")

            patt = Chem.MolFromSmarts('[C;!R]C=[C;!R]')
            if patt is None:
                raise ValueError("无效的反应模式")

            repl = Chem.MolFromSmiles(ffkm)
            if repl is None:
                raise ValueError("无效的FFKM分子结构")

            rms = AllChem.ReplaceSubstructs(
                mol, patt, repl,
                replaceAll=True,
                replacementConnectionPoint=connection_point
            )

            if not rms:
                raise ValueError("交联反应失败")

            return Chem.MolToSmiles(rms[0])

        except Exception as e:
            raise ValueError(f"交联反应模拟失败: {str(e)}")

    @staticmethod
    def validate_molecule(smiles: str) -> bool:
        """验证分子结构的有效性。

        Args:
            smiles: 分子的SMILES字符串

        Returns:
            bool: 如果分子结构有效返回True，否则返回False
        """
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    @staticmethod
    def optimize_structure(smiles: str) -> str:
        """优化分子结构。

        Args:
            smiles: 分子的SMILES字符串

        Returns:
            str: 优化后的SMILES字符串

        Raises:
            ValueError: 如果分子结构无效或优化失败
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("无效的分子结构")

        try:
            # 使用MMFF94s力场进行结构优化
            AllChem.MMFFOptimizeMolecule(mol)
            return Chem.MolToSmiles(mol)
        except Exception as e:
            raise ValueError(f"结构优化失败: {str(e)}")

    @staticmethod
    def smile2xyz(xyz_name: str, smile: str, randomSeed: int = None) -> None:
        """将SMILES字符串转换为XYZ格式文件，并创建相应的目录结构。

        Args:
            xyz_name: xyz文件名（包含.xyz后缀）
            smile: SMILES字符串
            randomSeed: 3D构象生成的随机种子（可选）

        Example:
            smile2xyz('CH4_0.xyz','C')
            将创建CH4目录并在其中生成CH4_0.xyz文件

        Note:
            xyz文件名应以_n.xyz结尾，其中n为索引（应为1或2）

        Raises:
            ValueError: 如果分子结构无效或文件操作失败
        """
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            raise ValueError("无效的分子结构")

        try:
            mol = Chem.AddHs(mol)
            if randomSeed is not None:
                AllChem.EmbedMolecule(mol, randomSeed=randomSeed)
            else:
                AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)

            xyz_str = Chem.MolToXYZBlock(mol)
            dir_name = xyz_name.split('_')[0]
            
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
                
            file_path = os.path.join(dir_name, xyz_name)
            with open(file_path, "w") as file:
                file.write(xyz_str)

        except Exception as e:
            raise ValueError(f"XYZ文件生成失败: {str(e)}")