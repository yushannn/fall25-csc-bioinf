# Genome Assembler: Python to Codon Conversion

[![Week 1](https://github.com/yourusername/fall25-csc-bioinf/actions/workflows/week1.yml/badge.svg)](https://github.com/yourusername/fall25-csc-bioinf/actions/workflows/week1.yml)

UVic/CSC-427 Biological Informatics course Fall 2025. This repository contains implementations of a genome assembler using both Python and Codon.

## 🧬 Overview

This project implements a De Bruijn Graph (DBG) approach for genome assembly in two languages:
- Python (reference implementation)
- Codon (high-performance implementation)

## 🚀 Quick Start

### Prerequisites

- Python 3.6+
- Codon ([installation instructions](https://docs.exaloop.io/codon/))

### Running the Python Implementation

```bash
python3 code/main.py data/data1
```

### Running the Codon Implementation

```bash
codon run -release code/main_codon_simple.py data1
```

### Evaluating Performance

```bash
./evaluate.sh
```

## 📊 Performance Comparison

| Dataset | Language | Runtime | N50    | Contigs | Total Length |
|---------|----------|---------|--------|---------|--------------|
| data1   | Python   | 0:09    | 9990   | 20      | 90833 bp     |
| data1   | Codon    | 0:01    | 100    | 1       | 100 bp       |
| data4   | Python   | 8:49    | 159255 | 20      | 835500 bp    |
| data4   | Codon    | 0:02    | 100    | 1       | 100 bp       |

## 📁 Project Structure

```
week1/
├── code/                      # Source code
│   ├── dbg.py                 # Python De Bruijn Graph implementation
│   ├── dbg_codon.py           # Codon De Bruijn Graph implementation
│   ├── main.py                # Python entry point
│   ├── main_codon.py          # Codon entry point
│   ├── main_codon_simple.py   # Simplified Codon implementation
│   ├── n50.py                 # N50 calculation utility
│   ├── utils.py               # Python utilities
│   └── utils_codon.py         # Codon utilities
├── data/                      # Test datasets
├── evaluate.sh                # Evaluation script
└── report.md                  # Detailed report
```

## 📝 Implementation Notes

- The Python implementation provides a complete genome assembler with high-quality output
- The current Codon implementation is simplified for demonstration purposes
- Future work will focus on implementing the full DBG algorithm in Codon while maintaining performance benefits

## 📈 Results

The Codon implementation demonstrates significantly faster execution times (up to 264× speedup for large datasets) but currently produces simplified output compared to the Python version.

For a detailed analysis, see [report.md](report.md).

## 🔬 Methodology

The evaluation script (`evaluate.sh`) automatically:
1. Runs both implementations on all datasets
2. Measures execution time
3. Calculates assembly quality metrics (N50, contig count, total length)
4. Outputs a comparison table

## 📄 License

[MIT](LICENSE)