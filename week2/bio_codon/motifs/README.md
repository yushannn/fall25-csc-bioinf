# Bio.motifs 模块详细文档

## 模块概述

这个目录包含Bio.motifs模块的核心实现，提供序列motif分析的完整功能。

## 📄 文件详情

### `__init__.py` (561行)
**功能**：主要的motif类和工具函数
- `Motif` 类：核心motif对象
- `SimpleAlignment` 类：序列对齐处理  
- `create()` 函数：从序列创建motif
- 常量定义和工具函数

**关键类和函数**：
```python
class Motif:
    """主要的motif类"""
    def __init__(self, alphabet=None, instances=None, counts=None)
    def consensus(self) -> str
    def anticonsensus(self) -> str  
    def reverse_complement(self) -> 'Motif'
    def format(self, format_type: str) -> str

class SimpleAlignment:
    """简单序列对齐类"""
    def __init__(self, sequences: List[str])
    def __getitem__(self, index: int) -> str
    def __len__(self) -> int

def create(instances, alphabet=None) -> Motif:
    """从序列列表创建motif对象"""
```

### `matrix.py` (547行)  
**功能**：各种矩阵类的实现
- `GenericPositionMatrix`：基础矩阵类
- `FrequencyPositionMatrix`：频率/计数矩阵
- `PositionWeightMatrix`：位置权重矩阵
- `PositionSpecificScoringMatrix`：位置特异性评分矩阵

**关键类**：
```python
class GenericPositionMatrix(dict):
    """基础位置矩阵类"""
    def __init__(self, alphabet: str, values: Dict[str, List[float]])
    @property
    def consensus(self) -> str
    @property  
    def anticonsensus(self) -> str

class FrequencyPositionMatrix(GenericPositionMatrix):
    """频率位置矩阵（计数矩阵）"""
    def normalize(self, pseudocounts=None) -> 'PositionWeightMatrix'
    def reverse_complement(self) -> 'FrequencyPositionMatrix'

class PositionWeightMatrix(GenericPositionMatrix):
    """位置权重矩阵（概率矩阵）"""
    def log_odds(self, background=None) -> 'PositionSpecificScoringMatrix' 
    def reverse_complement(self) -> 'PositionWeightMatrix'

class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """位置特异性评分矩阵"""
    def calculate(self, sequence: str) -> List[float]
    def search(self, sequence: str, threshold: float = 0.0) -> List[Tuple[int, float]]
    def reverse_complement(self) -> 'PositionSpecificScoringMatrix'
```

### `minimal.py` (456行)
**功能**：MEME minimal格式支持
- MEME格式文件的读写
- motif序列化和反序列化
- 与标准MEME工具的兼容性

**关键函数**：
```python
def read(handle) -> List[Motif]:
    """从MEME格式文件读取motifs"""

def write(motifs: List[Motif]) -> str:
    """将motifs写入MEME格式字符串"""

class Record:
    """MEME记录类"""
    def __init__(self)
    @property
    def motifs(self) -> List[Motif]
```

### `thresholds.py` (203行)
**功能**：阈值计算和分析
- 分数分布计算
- 阈值优化
- 统计分析工具

**关键函数**：
```python
def score_distribution(pssm, background=None, precision=10**4) -> Dict[float, float]:
    """计算PSSM分数分布"""

def threshold_from_p_value(pssm, p_value: float, background=None) -> float:
    """根据p值计算阈值"""

def p_value_from_threshold(pssm, threshold: float, background=None) -> float:
    """根据阈值计算p值"""
```

## 🔗 模块间关系

```
__init__.py
    ├── 导入并使用 matrix.py 中的矩阵类
    ├── 使用 minimal.py 进行格式化输出
    └── 使用 thresholds.py 进行阈值计算
    
matrix.py
    ├── 被 __init__.py 中的 Motif 类使用
    └── 为其他模块提供核心矩阵功能
    
minimal.py
    ├── 使用 __init__.py 中的 Motif 类
    └── 使用 matrix.py 中的矩阵类
    
thresholds.py
    └── 使用 matrix.py 中的 PSSM 类
```

## 🎯 设计原则

1. **API兼容性**：保持与BioPython的API兼容
2. **类型安全**：全面使用类型注解
3. **错误处理**：详细的参数验证和错误消息
4. **性能优化**：避免不必要的计算和内存分配
5. **可测试性**：每个功能都有对应的测试用例

## 🔧 扩展指南

要添加新功能：

1. **新的矩阵类型**：继承 `GenericPositionMatrix`
2. **新的文件格式**：参考 `minimal.py` 的实现
3. **新的分析工具**：添加到相应的模块中
4. **保持兼容性**：确保不破坏现有API

## 📊 代码统计

| 文件 | 行数 | 类数 | 函数数 | 主要功能 |
|------|------|------|--------|----------|
| `__init__.py` | 561 | 2 | 12+ | 核心motif功能 |
| `matrix.py` | 547 | 4 | 20+ | 矩阵计算 |
| `minimal.py` | 456 | 2 | 8+ | MEME格式 |
| `thresholds.py` | 203 | 0 | 6+ | 阈值分析 |
| **总计** | **1767** | **8** | **46+** | **完整motif工具集** |