# Week 2: Bio.motifs Codon Port

A Codon-compatible port of BioPython's Bio.motifs package for sequence motif analysis.

## Quick Start

```bash
# Clone repository
git clone https://github.com/yushannn/fall25-csc-bioinf.git
cd fall25-csc-bioinf/week2

# Run tests
python test.py    # All 17 tests passed
codon test.py     # All 17 tests passed
```

## Basic Usage

```python
from bio_codon import motifs

# Create motif from sequences
sequences = ["ACGT", "ACGG", "ACGA", "ACGC"]
motif = motifs.create(sequences)

# Access matrices
print(f"Consensus: {motif.consensus}")  # ACGN
fpm = motif.counts       # Frequency matrix
pwm = motif.pwm          # Weight matrix
pssm = motif.pssm        # Scoring matrix

# Score and search
score = pssm.calculate("ACGT")
matches = pssm.search("ATCGACGTACGA", threshold=0.0)
```

## Implementation Status

| Module | File | Lines | Status |
|--------|------|-------|--------|
| Core Motifs | `__init__.py` | 561 | ✓ Complete |
| Matrices | `matrix.py` | 547 | ✓ Complete |
| MEME Parser | `minimal.py` | 456 | ✓ Complete |
| Thresholds | `thresholds.py` | 203 | ✓ Complete |
| **Total** | | **1767** | |

## Features

**Implemented:**
- Motif creation and manipulation
- Frequency/Weight/Scoring matrices
- MEME format reading and writing
- Sequence scoring and searching
- Reverse complement operations
- Threshold calculations

**Not Implemented:**
- JASPAR, TRANSFAC formats
- Advanced statistical methods
- Motif visualization
- NumPy integration

## Testing

**17 tests** covering:
- Basic operations (5 tests)
- Extended functionality (3 tests)
- Edge cases (4 tests)
- Error handling (3 tests)
- Integration (2 tests)

All tests pass in both Python and Codon environments.

## Project Structure

```
week2/
├── README.md           # This file
├── report.md           # Detailed implementation report
├── ai.md              # AI usage documentation
├── test.py            # Test suite (17 tests)
└── bio_codon/         # Source code (1767 lines)
    └── motifs/
        ├── __init__.py        # Core classes
        ├── matrix.py          # Matrix operations
        ├── minimal.py         # MEME parser
        └── thresholds.py      # Threshold calculations
```

## Example Workflows

### Motif Analysis

```python
# Create and analyze motif
sequences = ["TACAA", "TACGC", "TACAC", "TACCC"]
motif = motifs.create(sequences)

print(f"Length: {motif.length}")
print(f"Consensus: {motif.consensus}")

# Get different matrices
counts = motif.counts
pwm = counts.normalize()
pssm = pwm.log_odds()
```

### Sequence Search

```python
# Search for matches
target = "ATGTACAATACGCATACAC"
matches = motif.pssm.search(target, threshold=0.0)

for start, end, score in matches:
    print(f"Match at {start}-{end}: score {score:.2f}")
```

### Reverse Complement

```python
# Work with both strands
forward = motifs.create(["ATCG", "ATGG"])
reverse = forward.reverse_complement()

print(f"Forward: {forward.consensus}")
print(f"Reverse: {reverse.consensus}")
```

### MEME Format

```python
from bio_codon.motifs import minimal

# Read MEME file
with open('motifs.meme', 'r') as f:
    record = minimal.read(f)

# Write MEME format
output = minimal.write([motif])
```

## Technical Details

**Codon Compatibility:**
- No NumPy dependencies
- Native Python types (dict, list)
- Full type annotations
- Cross-platform testing

**API Compatibility:**
- Matches BioPython function signatures
- Same property names and return types
- Compatible class hierarchies

**Error Handling:**
- Type validation at boundaries
- Descriptive error messages
- Range checking for all inputs

## Requirements

- Python 3.8+ (for Python mode)
- BioPython (for testing only)
- Codon v0.16+ (optional, for Codon mode)

## Documentation

- **report.md** - Complete technical documentation with implementation details, design decisions, and usage examples
- **ai.md** - Documentation of AI tool usage during development
- **test.py** - Test suite with 17 tests

## Known Limitations

1. Only MEME minimal format supported
2. Simplified threshold calculations
3. Not benchmarked for performance
4. Subset of BioPython features
5. Basic test coverage (17 tests)

## Development Info

- **Time:** 18-20 hours
- **Tests:** 17 (100% pass rate)
- **Code:** 1767 lines across 4 modules
- **AI Tools:** GitHub Copilot, ChatGPT (see ai.md)

## Notes

This is a course assignment implementing core motif analysis functionality. It is suitable for educational purposes and basic analysis, but not intended as a production replacement for BioPython.

For detailed technical information, see **report.md**.