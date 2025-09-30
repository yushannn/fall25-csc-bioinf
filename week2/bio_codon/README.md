# Bio.motifs Codon移植版本

这是Bio.motifs模块的Codon移植版本，实现了序列motif分析的核心功能。

## 🎯 移植目标

将以下BioPython模块移植到Codon：
- `Bio.motifs.__init__.py` → `bio_codon/motifs/__init__.py` ✅
- `Bio.motifs.matrix` → `bio_codon/motifs/matrix.py` ✅  
- `Bio.motifs.minimal` → `bio_codon/motifs/minimal.py` ✅
- `Bio.motifs.thresholds` → `bio_codon/motifs/thresholds.py` ✅

## 📁 目录结构

```
bio_codon/
├── __init__.py                 # 包初始化文件
├── README.md                   # 本文件
└── motifs/                     # motifs子模块
    ├── __init__.py            # 主motif类和工具函数
    ├── matrix.py              # 矩阵类 (FPM/PWM/PSSM)
    ├── minimal.py             # MEME格式支持
    ├── thresholds.py          # 阈值计算工具
    └── README.md              # 详细模块文档
```

## ✨ 已实现功能

### 核心类
- ✅ `Motif` - 主要的motif类
- ✅ `SimpleAlignment` - 序列对齐类  
- ✅ `create()` - motif创建函数

### 矩阵类
- ✅ `GenericPositionMatrix` - 基础矩阵类
- ✅ `FrequencyPositionMatrix` - 频率矩阵
- ✅ `PositionWeightMatrix` - 权重矩阵
- ✅ `PositionSpecificScoringMatrix` - 评分矩阵

### 格式支持
- ✅ MEME minimal格式读写
- ✅ 基本的motif格式化

### 工具功能
- ✅ 反向互补序列处理
- ✅ 阈值计算
- ✅ 序列搜索和评分

## 🚀 使用示例

```python
# 导入Codon版本
from bio_codon import motifs

# 创建motif
sequences = ["ACGT", "ACGG", "ACGA", "ACGC"]
m = motifs.create(sequences)

# 计算频率矩阵
fpm = m.counts.normalize()

# 计算权重矩阵  
pwm = fpm.log_odds()

# 搜索序列
scores = pwm.calculate("ACGTACGT")
```

## 🔄 与BioPython兼容性

| 功能 | BioPython | Codon移植 | 状态 |
|------|-----------|-----------|------|
| 基本motif创建 | ✅ | ✅ | 完全兼容 |
| 矩阵计算 | ✅ | ✅ | 完全兼容 |
| MEME格式 | ✅ | ✅ | 基本兼容 |
| 阈值计算 | ✅ | ⚠️ | 部分兼容* |
| NumPy集成 | ✅ | ❌ | 不支持** |

*部分BioPython的高级阈值功能未实现  
**Codon不支持NumPy，使用原生数学运算替代

## 🧪 测试

运行测试：
```bash
# Python环境
python test.py

# Codon环境  
codon test.py
```

## 📋 已知限制

1. **NumPy依赖**：移除了所有NumPy依赖，使用纯Codon实现
2. **高级格式**：只支持基本的MEME格式，不支持复杂的生物信息学格式
3. **性能**：某些矩阵运算可能比NumPy版本慢
4. **外部工具**：不支持调用外部生物信息学工具

## 🔧 开发说明

- 保持与BioPython API的兼容性
- 使用纯Codon代码，避免Python特定库
- 完整的类型注解支持
- 详细的错误处理和参数验证