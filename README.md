# Electret Calculator

一个用于计算电解质分子性质的Python工具包，专注于分子离子势能（IP）的计算。

## 功能特性

- 支持从SMILES格式输入分子结构
- 自动化的分子IP计算流程
- 支持正电荷和中性电荷状态的IP计算
- 集成Gaussian计算引擎
- 提供完整的数据处理和可视化工具

## 安装说明

### 依赖要求

- Python 3.7+
- RDKit
- Pandas
- Gaussian 09/16

### 安装步骤

```bash
git clone https://github.com/yourusername/electret_cal.git
cd electret_cal
pip install -r requirements.txt
```

## 使用方法

### 计算流程

完整的计算流程如下：

```
SMILES -> XYZ₀ -> GJF₀ -> SH₀ -> LOG₀ -> XYZ₁ -> GJF₁ -> SH₁ -> LOG₁ -> IP
```

其中：
- XYZ：分子的3D结构文件
- GJF：Gaussian输入文件
- SH：计算脚本
- LOG：计算结果日志

### 核心模块

项目包含以下核心模块：

- `smiles2xyz`：将SMILES转换为XYZ格式
- `xyz2gjf`：生成Gaussian输入文件
- `gjf2sh`：创建计算脚本
- `log2xyz`：提取优化后的结构
- `log2energy`：提取能量数据

### 计算类型

系统支持两种计算模式：

1. **正电荷状态计算** (`flow_pos`)
   - 计算分子在带正电荷状态下的IP
   - 适用于研究电子得失过程

2. **中性状态计算** (`flow_neu`)
   - 计算分子在中性状态下的IP
   - 用于基态性质研究

## API文档

### 主要函数

```python
flow_pos(smiles: str) -> float
```
计算分子在正电荷状态下的IP值

```python
flow_neu(smiles: str) -> float
```
计算分子在中性状态下的IP值

## 贡献指南

欢迎提交问题和改进建议！请遵循以下步骤：

1. Fork本仓库
2. 创建您的特性分支
3. 提交您的改动
4. 推送到您的分支
5. 创建Pull Request

## 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 联系方式

如有问题或建议，请通过以下方式联系我们：

- 提交Issue
- 发送邮件至：[yuhangu59@gmail.com]
