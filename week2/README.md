# Week 2: Bio.motifs Codon Port

A Codon-compatible implementation of BioPython's Bio.motifs package for sequence motif analysis.

## Project Structure

```
week2/
├── README.md                   # This file
├── test.py                     # Test suite (17 tests)
├── report.md                   # Implementation report
├── ai.md                       # AI usage documentation
└── bio_codon/                  # Codon port implementation
    ├── __init__.py            # Package initialization
    └── motifs/                # Motifs module
        ├── __init__.py        # Core Motif class (561 lines)
        ├── matrix.py          # Matrix classes (547 lines)
        ├── minimal.py         # MEME parser (456 lines)
        └── thresholds.py      # Threshold calculations (203 lines)
```

## Implementation Status

| BioPython Module | Codon Port | Status | Lines |
|------------------|------------|--------|-------|
| `Bio.motifs.__init__.py` | `bio_codon/motifs/__init__.py` | ✓ Complete | 561 |
| `Bio.motifs.matrix` | `bio_codon/motifs/matrix.py` | ✓ Complete | 547 |
| `Bio.motifs.minimal` | `bio_codon/motifs/minimal.py` | ✓ Complete | 456 |
| `Bio.motifs.thresholds` | `bio_codon/motifs/thresholds.py` | ✓ Complete | 203 |
| **Total** | | | **1767** |

## Quick Start

### Installation & Testing

```bash
# Clone repository
git clone https://github.com/yushannn/fall25-csc-bioinf.git
cd fall25-csc-bioinf/week2

# Run tests with Python
python test.py
# Output: All 17 tests passed

# Run tests with Codon (if installed)
codon test.py
# Output: All 17 tests passed
```

### Basic Usage

```python
from bio_codon import motifs

# Create a motif from sequences
sequences = ["ACGT", "ACGG", "ACGA", "ACGC"]
motif = motifs.create(sequences)

# Access motif properties
print(f"Length: {motif.length}")        # 4
print(f"Consensus: {motif.consensus}")  # ACGN

# Get different matrix representations
freq_matrix = motif.counts              # Frequency Position Matrix
weight_matrix = freq_matrix.normalize() # Position Weight Matrix
scoring_matrix = weight_matrix.log_odds() # Position Specific Scoring Matrix

# Score a sequence
score = scoring_matrix.calculate("ACGT")
print(f"Score: {score}")
```

## Features

### Implemented ✓

**Core Functionality:**
- Motif creation from sequence lists
- Frequency Position Matrix (FPM)
- Position Weight Matrix (PWM)
- Position Specific Scoring Matrix (PSSM)
- Consensus sequence calculation
- Reverse complement (DNA/RNA)

**Matrix Operations:**
- Normalization with pseudocounts
- Log-odds calculations
- Sequence scoring
- Motif searching with thresholds

**File Format Support:**
- MEME minimal format reading
- MEME minimal format writing

**Statistical Analysis:**
- Score distribution calculation
- Threshold determination (FPR/FNR)
- P-value calculations

### Not Implemented ✗

- JASPAR format parser
- TRANSFAC format parser
- MAST format parser
- AlignACE format parser
- Motif comparison functions
- Motif visualization
- NumPy integration (Codon limitation)

## Testing

**Test Suite: 17 tests**

| Category | Count | Coverage |
|----------|-------|----------|
| Basic Operations | 5 | Motif creation, matrix operations, consensus |
| Extended Functionality | 3 | Reverse complement, PSSM scoring, search |
| Edge Cases | 4 | Single/short/long sequences, large datasets |
| Error Handling | 3 | Empty input, mixed lengths, type errors |
| Integration | 2 | Full workflow, invalid characters |

**Running Tests:**
```bash
# Python
python test.py

# Codon
codon test.py
```

Both should output: `All 17 tests passed`

## API Documentation

### Creating Motifs

```python
from bio_codon import motifs

# From sequence list
sequences = ["TACAA", "TACGC", "TACAC", "TACCC"]
motif = motifs.create(sequences)

# From counts dictionary
counts = {
    'A': [4, 0, 4, 0, 4],
    'C': [0, 4, 0, 4, 0],
    'G': [0, 0, 0, 0, 0],
    'T': [0, 0, 0, 0, 0]
}
motif = motifs.Motif(alphabet="ACGT", counts=counts)
```

### Matrix Operations

```python
# Get frequency matrix
fpm = motif.counts

# Normalize to get PWM
pwm = fpm.normalize(pseudocounts={'A': 0.5, 'C': 0.5, 'G': 0.5, 'T': 0.5})

# Calculate log-odds for PSSM
background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
pssm = pwm.log_odds(background=background)

# Or use shortcuts
pwm = motif.pwm
pssm = motif.pssm
```

### Sequence Scoring

```python
# Score a single position
score = pssm.calculate("TACAA", start=0)

# Search entire sequence
matches = pssm.search("ATACTACAATACGG", threshold=5.0)
for start, end, score in matches:
    print(f"Position {start}-{end}: score {score:.2f}")
```

### Reverse Complement

```python
# DNA motif
dna_motif = motifs.create(["ATCG", "ATGG"])
rc_motif = dna_motif.reverse_complement()

print(dna_motif.consensus)  # Original
print(rc_motif.consensus)   # Reverse complement
```

### MEME Format

```python
from bio_codon.motifs import minimal

# Read MEME file
with open('motifs.meme', 'r') as f:
    record = minimal.read(f)

# Access motifs
for motif in record:
    print(f"Motif: {motif.name}, Length: {motif.length}")

# Write MEME format
output = minimal.write([motif])
print(output)
```

## Technical Details

### Codon Compatibility

**Adaptations Made:**
- Removed NumPy dependencies (Codon doesn't support NumPy)
- Used native Python types (dict, list) for matrices
- Adapted exception handling for Codon constraints
- Added comprehensive type annotations

**Type Annotations:**
```python
def create(instances: List[str], alphabet: str = "ACGT") -> Motif
def normalize(self, pseudocounts: Optional[Dict[str, float]] = None) -> PositionWeightMatrix
def calculate(self, sequence: str, start: int = 0) -> float
```

### Numerical Stability

**Handling Edge Cases:**
```python
MIN_PROBABILITY = 1e-10  # Avoid log(0)

# In log-odds calculation
frequency = max(self[letter][i], MIN_PROBABILITY)
bg_freq = max(background[letter], MIN_PROBABILITY)
log_odds = math.log2(frequency / bg_freq)
```

### Error Handling

**Input Validation:**
- Type checking (list vs string vs int)
- Length consistency (all sequences same length)
- Alphabet validation (only valid characters)
- Range checking (non-negative values)

**Example:**
```python
# This will raise ValueError
motifs.create(["ATCG", "AT"])  # Different lengths

# This will raise TypeError  
motifs.create("ATCG")  # Not a list

# This will raise ValueError
motifs.create(["ATXG"])  # Invalid character X
```

## Functionality Comparison

| Feature | BioPython | This Port | Compatibility |
|---------|-----------|-----------|---------------|
| Motif creation | ✓ | ✓ | 100% |
| FPM, PWM, PSSM | ✓ | ✓ | 100% |
| Consensus sequence | ✓ | ✓ | 100% |
| Reverse complement | ✓ | ✓ | 100% |
| Sequence scoring | ✓ | ✓ | 100% |
| Motif search | ✓ | ✓ | 100% |
| MEME format | ✓ | ✓ | Minimal only |
| JASPAR format | ✓ | ✗ | Not implemented |
| TRANSFAC format | ✓ | ✗ | Not implemented |
| Threshold calculation | ✓ | ✓ | Simplified |
| NumPy integration | ✓ | ✗ | Codon limitation |

## Known Limitations

1. **Format Support**
   - Only MEME minimal format implemented
   - JASPAR, TRANSFAC, MAST not supported

2. **Statistical Methods**
   - Threshold calculations are simplified
   - Some advanced statistics not implemented

3. **Performance**
   - Not benchmarked against BioPython
   - No performance optimization done yet

4. **Features**
   - Subset of full BioPython.motifs
   - Some advanced methods missing

5. **Testing**
   - 17 tests cover basic functionality
   - Not exhaustive test coverage

## Example Workflows

### Complete Analysis Pipeline

```python
from bio_codon import motifs

# 1. Create motif
sequences = ["TACAA", "TACGC", "TACAC", "TACCC"]
motif = motifs.create(sequences)

# 2. Examine properties
print(f"Motif length: {motif.length}")
print(f"Consensus: {motif.consensus}")

# 3. Get matrices
print("Frequency matrix:")
print(motif.counts)

print("\nPosition Weight Matrix:")
print(motif.pwm)

print("\nPosition Specific Scoring Matrix:")
print(motif.pssm)

# 4. Score sequences
test_sequences = ["TACAA", "TACGC", "ATGCG"]
for seq in test_sequences:
    score = motif.pssm.calculate(seq)
    print(f"{seq}: {score:.2f}")

# 5. Search for matches
target = "ATGTACAATACGCATACAC"
matches = motif.pssm.search(target, threshold=0.0)
print(f"\nFound {len(matches)} matches in target sequence")
```

### Working with Thresholds

```python
from bio_codon.motifs import thresholds

# Create score distribution
distribution = thresholds.score_distribution(motif)

# Calculate threshold for 5% false positive rate
threshold = distribution.threshold_fpr(0.05)
print(f"Threshold for 5% FPR: {threshold:.2f}")

# Calculate FPR for a given score
fpr = distribution.score_to_fpr(5.0)
print(f"FPR at score 5.0: {fpr:.4f}")

# Find balanced threshold
balanced = distribution.threshold_balanced()
print(f"Balanced threshold: {balanced:.2f}")
```

### Strand-Specific Analysis

```python
# Create motif
forward_motif = motifs.create(["ATCGATCG", "ATCGATCG"])

# Get reverse complement
reverse_motif = forward_motif.reverse_complement()

# Search both strands
sequence = "ATCGATCGATCGATCG"

forward_matches = forward_motif.pssm.search(sequence)
reverse_matches = reverse_motif.pssm.search(sequence)

print(f"Forward strand: {len(forward_matches)} matches")
print(f"Reverse strand: {len(reverse_matches)} matches")
```

## Development Information

### Time Investment
- **Total:** 18-20 hours
- Code structure: 5 hours
- Matrix operations: 4 hours
- Format parsing: 4 hours
- Testing: 4 hours
- Documentation: 1 hour

### Tools Used
- GitHub Copilot for code suggestions
- ChatGPT for understanding BioPython APIs
- See `ai.md` for detailed AI usage documentation

### Design Philosophy
1. **Simplicity:** Use native Python types
2. **Compatibility:** Match BioPython API
3. **Safety:** Comprehensive input validation
4. **Clarity:** Extensive type annotations and docstrings

## Requirements

**For Python mode:**
- Python 3.8+
- BioPython (`pip install biopython`)

**For Codon mode:**
- Codon compiler v0.16+
- No additional dependencies

**For CI:**
- See `.github/workflows/week2.yml`

## Files

- `test.py` - 17 tests for both Python and Codon
- `report.md` - Detailed implementation report
- `ai.md` - AI assistance documentation
- `bio_codon/` - Source code (1767 lines)

## Future Enhancements

**Potential additions (not implemented):**
- JASPAR format parser
- TRANSFAC format parser
- Motif comparison algorithms
- Performance benchmarking
- Additional statistical methods
- Extended test coverage
- Visualization (if graphics support added)

## CI/CD

GitHub Actions workflow runs:
- Python 3.13 tests
- Codon installation (optional)
- All 17 tests must pass

See `.github/workflows/week2.yml` for configuration.

## Notes

This is a course assignment for porting BioPython functionality to Codon. It implements core motif analysis features but is not intended as a complete replacement for BioPython. Use BioPython for production work requiring full feature support.

## License

This is educational code created as a course assignment. Original BioPython code is under the Biopython License Agreement.

## Acknowledgments

- BioPython developers for the original implementation
- Course instructor for assignment design
- GitHub Copilot and ChatGPT for development assistance