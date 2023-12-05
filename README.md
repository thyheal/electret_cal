# electret_cal
This is a repo store script to calculate some properties of electret
# 分子IP计算工具文档简介

这是一些用于计算分子离子势能（IP）的工具。计算流程可以简要概括如下：

1. **计算流程**：

smiles -> xyz0 -> gjf0 -> sh0 -> log0 -> xyz1 -> gjf1 -> sh1 -> log1 -> ip


在计算过程中，使用了一些模块来生成所需的文件，如`smiles2xyz`、`xyz2gjf`、`gjf2sh`、`log2xyz`和`log2energy`。

2. **计算过程分为两部分**：

- 第一部分：计算正电荷状态下的IP。
- 第二部分：计算中性电荷状态下的IP。

因此，流程被划分为两个函数：`flow_pos` 用于正电荷状态的计算，`flow_neu` 用于中性电荷状态的计算。

这个文档简介提供了有关分子IP计算工具的高级概述，以及计算流程的基本说明。详细的工具使用和计算方法可以在后续文档或说明中找到。
