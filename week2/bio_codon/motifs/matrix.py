"""Support for various forms of sequence motif matrices (Codon port).

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""

from typing import Dict, List, Optional, Union, Any
import math


class GenericPositionMatrix(dict):
    """Base class for the support of position matrix operations."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        self.length = None
        for letter in alphabet:
            if self.length is None:
                self.length = len(values[letter])
            elif self.length != len(values[letter]):
                raise Exception("data has inconsistent lengths")
            # Convert values to Python floats
            self[letter] = [float(v) for v in values[letter]]
        self.alphabet = alphabet

    def __str__(self) -> str:
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [f"{value:6.2f}" for value in self[letter]]
            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text

    def __repr__(self) -> str:
        return str(self)

    def transpose(self) -> List[Dict[str, float]]:
        """Return the transpose of the matrix as a list of dicts."""
        transposed = []
        for i in range(self.length):
            position = {}
            for letter in self.alphabet:
                position[letter] = self[letter][i]
            transposed.append(position)
        return transposed


class FrequencyPositionMatrix(GenericPositionMatrix):
    """A frequency (count) matrix."""

    def __init__(self, alphabet: str, counts: Dict[str, List[int]]):
        """Initialize the FrequencyPositionMatrix."""
        # Convert counts to floats for compatibility
        float_counts = {}
        for letter in alphabet:
            if letter in counts:
                float_counts[letter] = [float(c) for c in counts[letter]]
            else:
                # If letter not in counts, initialize with zeros
                length = len(next(iter(counts.values()))) if counts else 0
                float_counts[letter] = [0.0] * length
        
        super().__init__(alphabet, float_counts)

    @property
    def counts(self) -> Dict[str, List[int]]:
        """Return the counts as integers."""
        return {letter: [int(c) for c in self[letter]] for letter in self.alphabet}

    def normalize(self, pseudocounts: Optional[Dict[str, float]] = None) -> "PositionWeightMatrix":
        """Normalize the frequency matrix to create a position weight matrix."""
        if pseudocounts is None:
            pseudocounts = {letter: 0.0 for letter in self.alphabet}
        
        # Calculate position weight matrix
        pwm_values = {}
        for letter in self.alphabet:
            pwm_values[letter] = []
            for i in range(self.length):
                # Add pseudocount and normalize
                total = sum(self[l][i] + pseudocounts[l] for l in self.alphabet)
                if total > 0:
                    frequency = (self[letter][i] + pseudocounts[letter]) / total
                else:
                    frequency = 1.0 / len(self.alphabet)  # Uniform distribution
                pwm_values[letter].append(frequency)
        
        return PositionWeightMatrix(self.alphabet, pwm_values)


class PositionWeightMatrix(GenericPositionMatrix):
    """A position weight matrix (PWM)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the PositionWeightMatrix."""
        super().__init__(alphabet, values)
        
        # Verify that each position sums to 1 (approximately)
        for i in range(self.length):
            total = sum(self[letter][i] for letter in self.alphabet)
            if abs(total - 1.0) > 1e-6:
                # Normalize if not already normalized
                for letter in self.alphabet:
                    if total > 0:
                        self[letter][i] /= total

    def log_odds(self, background: Optional[Dict[str, float]] = None) -> "PositionSpecificScoringMatrix":
        """Calculate log-odds scores to create a PSSM."""
        if background is None:
            # Use uniform background
            background = {letter: 1.0 / len(self.alphabet) for letter in self.alphabet}
        
        pssm_values = {}
        for letter in self.alphabet:
            pssm_values[letter] = []
            for i in range(self.length):
                frequency = self[letter][i]
                bg_freq = background[letter]
                if frequency > 0 and bg_freq > 0:
                    log_odds = math.log2(frequency / bg_freq)
                else:
                    log_odds = float('-inf')  # Very low score
                pssm_values[letter].append(log_odds)
        
        return PositionSpecificScoringMatrix(self.alphabet, pssm_values)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """A position-specific scoring matrix (PSSM)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the PositionSpecificScoringMatrix."""
        super().__init__(alphabet, values)

    @property
    def min(self) -> float:
        """Return the minimum possible score."""
        min_score = 0.0
        for i in range(self.length):
            position_min = min(self[letter][i] for letter in self.alphabet)
            min_score += position_min
        return min_score

    @property
    def max(self) -> float:
        """Return the maximum possible score."""
        max_score = 0.0
        for i in range(self.length):
            position_max = max(self[letter][i] for letter in self.alphabet)
            max_score += position_max
        return max_score

    def calculate(self, sequence: str, start: int = 0) -> float:
        """Calculate the score for a sequence at the given start position."""
        if start + self.length > len(sequence):
            raise ValueError("Sequence too short for motif")
        
        score = 0.0
        for i in range(self.length):
            nucleotide = sequence[start + i].upper()
            if nucleotide in self.alphabet:
                score += self[nucleotide][i]
            else:
                # Handle ambiguous nucleotides by using minimum score
                score += min(self[letter][i] for letter in self.alphabet)
        
        return score

    def search(self, sequence: str, threshold: float = 0.0) -> List[tuple]:
        """Search for motif matches in a sequence above the threshold."""
        matches = []
        seq_len = len(sequence)
        
        for start in range(seq_len - self.length + 1):
            score = self.calculate(sequence, start)
            if score >= threshold:
                end = start + self.length
                matches.append((start, end, score))
        
        return matches

    def reverse_complement(self) -> "PositionSpecificScoringMatrix":
        """Return the reverse complement of the PSSM."""
        # This assumes DNA alphabet
        complement_map = {"A": "T", "C": "G", "G": "C", "T": "A"}
        
        rc_values = {}
        for letter in self.alphabet:
            if letter in complement_map:
                complement_letter = complement_map[letter]
                rc_values[complement_letter] = list(reversed(self[letter]))
            else:
                # For non-standard letters, just reverse
                rc_values[letter] = list(reversed(self[letter]))
        
        return PositionSpecificScoringMatrix(self.alphabet, rc_values)

    def mean(self, background: Optional[Dict[str, float]] = None) -> float:
        """Calculate the mean score of the PSSM."""
        if background is None:
            background = {letter: 1.0 / len(self.alphabet) for letter in self.alphabet}
        
        mean_score = 0.0
        for i in range(self.length):
            position_mean = sum(background[letter] * self[letter][i] for letter in self.alphabet)
            mean_score += position_mean
        
        return mean_score