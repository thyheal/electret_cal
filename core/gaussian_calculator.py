from typing import Optional, Dict, Union
import os
from pathlib import Path

class GaussianCalculator:
    """高斯计算核心类，负责生成和执行高斯计算任务。

    该类处理高斯计算的配置、输入文件生成和计算执行。它是对原始gaussian_cal.py的重构版本，
    提供了更好的类型注解、错误处理和文档。

    Attributes:
        method (str): 计算方法
        basis (str): 基组
        charge (str): 电荷类型 ('pos', 'neg', 'neu')
        eps (float): PCM模型的介电常数
        opt (bool): 是否进行结构优化
        dispersion (bool): 是否包含色散校正
        polar (bool): 是否计算极化率
        volume (bool): 是否计算体积
        pcm (bool): 是否使用PCM溶剂化模型
        wfn (bool): 是否输出波函数文件
        debug (bool): 是否在调试模式下运行
    """

    def __init__(self,
                 method: str,
                 basis: str,
                 charge: str,
                 eps: float = 0,
                 opt: bool = False,
                 dispersion: bool = False,
                 polar: bool = False,
                 volume: bool = False,
                 pcm: bool = False,
                 wfn: bool = True,
                 debug: bool = False,
                 parent_dir: Optional[str] = None) -> None:
        """初始化GaussianCalculator实例。

        Args:
            method: 计算方法
            basis: 基组
            charge: 电荷类型 ('pos', 'neg', 'neu')
            eps: PCM模型的介电常数
            opt: 是否进行结构优化
            dispersion: 是否包含色散校正
            polar: 是否计算极化率
            volume: 是否计算体积
            pcm: 是否使用PCM溶剂化模型
            wfn: 是否输出波函数文件
            debug: 是否在调试模式下运行
        """
        self.method = method
        self.basis = basis
        self.charge = charge.lower()
        self.eps = eps
        self.opt = opt
        self.dispersion = dispersion
        self.polar = polar
        self.volume = volume
        self.pcm = pcm
        self.wfn = wfn
        self.debug = debug
        self.parent_dir = Path(parent_dir) if parent_dir else Path('.')

        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """验证输入参数的有效性。"""
        if self.charge not in ['pos', 'neg', 'neu']:
            raise ValueError(f"无效的电荷类型: {self.charge}. 必须是 'pos', 'neg' 或 'neu'")

    def run(self, xyz_name: str, dir_name: Optional[str] = None) -> None:
        """运行高斯计算任务。

        Args:
            xyz_name: XYZ文件名
            dir_name: 自定义目录名称（可选，如果不提供则使用xyz文件名的前缀）

        Raises:
            FileNotFoundError: 如果找不到必要的文件
            IOError: 如果文件操作失败
        """
        try:
            self._prepare_calculation(xyz_name, dir_name)
            self._generate_input_file()
            self._generate_shell_script()
            self._submit_job()
        except Exception as e:
            raise RuntimeError(f"计算任务执行失败: {str(e)}")

    def _prepare_calculation(self, xyz_name: str, dir_name: Optional[str] = None) -> None:
        """准备计算所需的文件和路径。

        Args:
            xyz_name: XYZ文件名
            dir_name: 自定义目录名称（可选，如果不提供则使用xyz文件名的前缀）
        """
        self.xyz_name = xyz_name
        self.dir_name = dir_name if dir_name is not None else Path(xyz_name).stem.split('_')[0]
        
        # 设置文件路径
        self.xyz_path = self.parent_dir / self.dir_name / self.xyz_name
        self.header_path = Path('.') / 'header' / 'gjf_header.txt'
        
        # 检查必要文件是否存在
        if not self.xyz_path.exists():
            raise FileNotFoundError(f"找不到XYZ文件: {self.xyz_path}")
        if not self.header_path.exists():
            raise FileNotFoundError(f"找不到头文件模板: {self.header_path}")

    def _generate_input_file(self) -> None:
        """生成高斯输入文件。"""
        charge_suffix = {'pos': 'pos', 'neg': 'neg', 'neu': 'neu'}[self.charge]
        self.gjf_path = self.parent_dir / self.dir_name / f"{self.dir_name}_{charge_suffix}.gjf"
        self.wfn_name = f"{self.dir_name}_{charge_suffix}.wfn"

        # 读取并修改模板
        with open(self.header_path, 'r') as f:
            header_content = f.read()

        # 替换模板中的占位符
        replacements = {
            '__title__': f"{self.dir_name}_{charge_suffix}",
            '__dir__': self.dir_name,
            '__method__': self.method,
            '__basis__': self.basis,
            '__charge__': {'pos': '1 2', 'neg': '-1 2', 'neu': '0 1'}[self.charge],
            '__opt__': 'opt' if self.opt else '',
            '__dispersion__': 'em=gd3' if self.dispersion else '',
            '__polar__': 'polar' if self.polar else '',
            '__Volume__': 'volume' if self.volume else '',
            '__PCM__': 'SCRF=(PCM,Solvent=Generic,Read)' if self.pcm else '',
            '__wfn__': 'out=wfn' if self.wfn else ''
        }

        for key, value in replacements.items():
            header_content = header_content.replace(key, value)

        # 读取原子坐标
        with open(self.xyz_path, 'r') as f:
            coordinates = f.readlines()[1:]

        # 写入完整的输入文件
        with open(self.gjf_path, 'w') as f:
            f.write(header_content)
            f.writelines(coordinates)
            f.write('\n')
            
            if self.pcm:
                f.write(f'eps={self.eps}\n\n')
            
            if self.wfn:
                f.write(f'{self.dir_name}/{self.wfn_name}\n')
            
            f.write('\n' * 4)

    def _generate_shell_script(self) -> None:
        """生成提交作业的Shell脚本。"""
        self.sh_path = self.gjf_path.with_suffix('.sh')
        sh_template_path = Path('.') / 'header' / 'sub_g09.txt'

        if not sh_template_path.exists():
            raise FileNotFoundError(f"找不到Shell脚本模板: {sh_template_path}")

        with open(sh_template_path, 'r') as f:
            sh_content = f.read()

        # 替换Shell脚本中的占位符
        replacements = {
            '__dir__': self.dir_name,
            '__name__': self.gjf_path.stem,
            '__changegjf__': str(self.gjf_path),
            '__changelog__': str(self.gjf_path.with_suffix('.log')),
            '__debug__': 'i8cpu' if self.debug else 'F1cpu'
        }

        for key, value in replacements.items():
            sh_content = sh_content.replace(key, value)

        with open(self.sh_path, 'w') as f:
            f.write(sh_content)

    def _submit_job(self) -> None:
        """提交计算作业到计算系统。"""
        os.system(f"sbatch {self.sh_path}")

    def generate_log(self) -> None:
        """生成计算日志文件。"""
        txt_path = self.parent_dir / 'A_intro.txt'
        header_path = Path('header') / 'A_intro.txt'

        if not header_path.exists():
            raise FileNotFoundError(f"找不到日志模板文件: {header_path}")

        with open(header_path, 'r') as f:
            content = f.read()

        # 替换日志文件中的占位符
        replacements = {
            '__method__': self.method,
            '__basis__': self.basis,
            '__charge__': {'pos': 'IP', 'neg': 'EA', 'neu': '0'}[self.charge],
            '__dispersion__': 'em=gd3' if self.dispersion else 'None',
            '__PCM__': 'SCRF=(PCM,Solvent=Generic,Read)' if self.pcm else 'None',
            '__EPS__': str(self.eps) if self.pcm else 'None',
            '__wfn__': 'out=wfn' if self.wfn else ''
        }

        for key, value in replacements.items():
            content = content.replace(key, str(value))

        with open(txt_path, 'w') as f:
            f.write(content)