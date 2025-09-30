# BioPython 源文件

这个目录包含从BioPython提取的原始源文件，用作Codon移植的参考。

## 文件结构

```
biopython-source/
├── README.md           # 本文件
├── __init__.py        # Bio.motifs.__init__.py 的原始版本
├── matrix.py          # Bio.motifs.matrix 的原始版本  
├── minimal.py         # Bio.motifs.minimal 的原始版本
└── thresholds.py      # Bio.motifs.thresholds 的原始版本
```

## 版本信息

- 来源：BioPython 1.84+
- 版权：Biopython Contributors
- 许可证：Biopython License Agreement 或 BSD 3-Clause License

## 用途

这些文件仅用作：

1. **参考**：了解BioPython的原始实现
2. **对比**：与Codon移植版本进行功能对比
3. **文档**：保留原始API文档和注释

## 重要说明

⚠️ **这些文件不会在Codon环境中运行**，因为它们依赖于：
- NumPy
- 其他BioPython模块
- Python特定的库

实际的Codon移植版本位于 `../bio_codon/` 目录中。