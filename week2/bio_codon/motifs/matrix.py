"""Support for various forms of sequence motif matrices (Codon port).

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""

from typing import Dict, List, Optional, Union, Any
import math

# Constants
DEFAULT_LOG_BASE = 2
MIN_PROBABILITY = 1e-10  # Minimum probability to avoid log(0)
UNIFORRM_BACKGROUND_DNA = 0.25  # Equal probability for A, C, G, T
DEFAULT_PRECISION = 6  # Number of decimal places for display
MAX_MATRIX_SIZE = 10000  # Safety limit for matrix dimensions

# DNA complement mapping for reverse complement
DNA_COMPLEMENT_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}


class GenericPositionMatrix(dict):
    """Base class for the support of position matrix operations.
    
    This class provides common functionality for position-specific matrices
    used in motif analysis, including frequency matrices, weight matrices,
    and scoring matrices.
    
    Attributes:
        alphabet: String containing valid characters (e.g., "ACGT").
        length: Number of positions in the matrix.
    
    The matrix data is stored as a dictionary where keys are alphabet
    characters and values are lists of numerical values for each position.
    """

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the position matrix.
        
        Args:
            alphabet: String of valid characters (e.g., "ACGT").
            values: Dictionary mapping each alphabet character to a list
                   of numerical values for each position. All lists must
                   have the same length.
        
        Raises:
            ValueError: If the value lists have inconsistent lengths, or if
                       required alphabet characters are missing from values.
            TypeError: If inputs have incorrect types.
        """
        if not isinstance(alphabet, str) or not alphabet:
            raise TypeError("alphabet must be a non-empty string")
        
        if not isinstance(values, dict):
            raise TypeError("values must be a dictionary")
        
        if not values:
            raise ValueError("values dictionary cannot be empty")
        
        # Check that all alphabet characters are present
        missing_chars = set(alphabet) - set(values.keys())
        if missing_chars:
            raise ValueError(f"Missing alphabet characters in values: {sorted(missing_chars)}")
        
        self.length = None
        for letter in alphabet:
            if letter not in values:
                raise ValueError(f"Alphabet character '{letter}' not found in values")
            
            if not isinstance(values[letter], (list, tuple)):
                raise TypeError(f"Values for '{letter}' must be a list or tuple")
            
            if self.length is None:
                self.length = len(values[letter])
                if self.length == 0:
                    raise ValueError("Value lists cannot be empty")
            elif self.length != len(values[letter]):
                raise ValueError(f"Inconsistent lengths: expected {self.length}, "
                               f"got {len(values[letter])} for character '{letter}'")
            
            # Convert values to Python floats with error handling
            try:
                self[letter] = [float(v) for v in values[letter]]
            except (ValueError, TypeError) as e:
                raise ValueError(f"Cannot convert values for '{letter}' to float: {e}")
        
        self.alphabet = alphabet

    def __str__(self) -> str:
        """Return a string containing nucleotides and values in the Matrix.
        
        Creates a formatted table showing position indices as column headers
        and alphabet characters as row labels, with matrix values displayed
        in a readable format.
        
        Returns:
            str: Multi-line string representation of the matrix.
        """
        # Use constants for formatting
        col_width = max(DEFAULT_PRECISION + 2, 6)  # Ensure minimum width
        
        # Header line with position indices
        words = [f"{i:{col_width}d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        
        # Data lines for each alphabet character
        for letter in self.alphabet:
            words = [f"{value:{col_width}.{DEFAULT_PRECISION-4}f}" for value in self[letter]]
            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        
        text = "\n".join(lines) + "\n"
        return text

    def __repr__(self) -> str:
        """Return string representation for debugging."""
        return str(self)

    def transpose(self) -> List[Dict[str, float]]:
        """Return the transpose of the matrix as a list of dicts.
        
        Converts the matrix from character-indexed format to position-indexed
        format, which is useful for position-specific analysis.
        
        Returns:
            List of dictionaries, where each dictionary represents one position
            and maps alphabet characters to their values at that position.
        
        Example:
            Matrix: {'A': [1.0, 2.0], 'C': [3.0, 4.0]}
            Returns: [{'A': 1.0, 'C': 3.0}, {'A': 2.0, 'C': 4.0}]
        """
        transposed = []
        for i in range(self.length):
            position = {}
            for letter in self.alphabet:
                position[letter] = self[letter][i]
            transposed.append(position)
        return transposed


class FrequencyPositionMatrix(GenericPositionMatrix):
    """A frequency (count) matrix.
    
    Represents the count of each alphabet character at each position
    in a collection of aligned sequences. This is the most basic form
    of motif representation, containing raw frequency data.
    
    The matrix stores integer counts internally as floating-point values
    for compatibility with mathematical operations, but provides access
    to integer counts through the counts property.
    """

    def __init__(self, alphabet: str, counts: Dict[str, List[int]]):
        """Initialize the FrequencyPositionMatrix.
        
        Args:
            alphabet: String of valid characters (e.g., "ACGT").
            counts: Dictionary mapping alphabet characters to lists of
                   integer counts. Missing characters are initialized with zeros.
        """
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
        """Return the counts as integers.
        
        Converts the internal floating-point representation back to
        integer counts for compatibility and display purposes.
        
        Returns:
            Dictionary mapping alphabet characters to lists of integer counts.
        """
        return {letter: [int(c) for c in self[letter]] for letter in self.alphabet}

    def normalize(self, pseudocounts: Optional[Dict[str, float]] = None) -> "PositionWeightMatrix":
        """Normalize the frequency matrix to create a position weight matrix.
        
        Converts raw counts to frequencies (probabilities) by dividing each
        count by the total count at that position. Pseudocounts can be added
        to avoid zero probabilities.
        
        Args:
            pseudocounts: Optional dictionary of pseudocounts to add to each
                         character before normalization. If None, uses zeros.
                         Helps avoid zero probabilities in sparse data.
        
        Returns:
            PositionWeightMatrix: Normalized matrix with frequencies summing to 1
                                 at each position.
        
        Raises:
            ValueError: If pseudocounts contain negative values or invalid characters.
            TypeError: If pseudocounts is not a dictionary.
        """
        if pseudocounts is None:
            pseudocounts = {letter: 0.0 for letter in self.alphabet}
        else:
            if not isinstance(pseudocounts, dict):
                raise TypeError("pseudocounts must be a dictionary")
            
            # Validate pseudocounts
            for letter in self.alphabet:
                if letter not in pseudocounts:
                    pseudocounts[letter] = 0.0
                elif not isinstance(pseudocounts[letter], (int, float)):
                    raise TypeError(f"Pseudocount for '{letter}' must be a number")
                elif pseudocounts[letter] < 0:
                    raise ValueError(f"Pseudocount for '{letter}' cannot be negative")
        
        # Calculate position weight matrix
        pwm_values = {}
        for letter in self.alphabet:
            pwm_values[letter] = []
            for i in range(self.length):
                # Add pseudocount and normalize
                total = sum(self[l][i] + pseudocounts.get(l, 0.0) for l in self.alphabet)
                if total > 0:
                    frequency = (self[letter][i] + pseudocounts.get(letter, 0.0)) / total
                else:
                    # Uniform distribution if all counts are zero
                    frequency = 1.0 / len(self.alphabet)
                
                pwm_values[letter].append(frequency)
        
        return PositionWeightMatrix(self.alphabet, pwm_values)


class PositionWeightMatrix(GenericPositionMatrix):
    """A position weight matrix (PWM).
    
    Represents the frequency (probability) of each alphabet character
    at each position in a motif. Each position sums to 1.0, making
    this suitable for probabilistic analysis and log-odds calculations.
    
    PWMs are created by normalizing frequency matrices and are used
    as an intermediate step in creating scoring matrices.
    """

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the PositionWeightMatrix.
        
        Args:
            alphabet: String of valid characters (e.g., "ACGT").
            values: Dictionary mapping alphabet characters to lists of
                   frequency values. Values at each position should sum to 1.0.
        
        Note:
            If values don't sum to 1.0 at each position, they will be
            automatically normalized during initialization.
        """
        super().__init__(alphabet, values)
        
        # Verify that each position sums to 1 (approximately)
        tolerance = MIN_PROBABILITY  # Use defined constant
        for i in range(self.length):
            total = sum(self[letter][i] for letter in self.alphabet)
            if abs(total - 1.0) > tolerance:
                # Normalize if not already normalized
                if total > MIN_PROBABILITY:  # Avoid division by very small numbers
                    for letter in self.alphabet:
                        self[letter][i] /= total
                else:
                    # If total is too small, use uniform distribution
                    uniform_prob = 1.0 / len(self.alphabet)
                    for letter in self.alphabet:
                        self[letter][i] = uniform_prob

    def log_odds(self, background: Optional[Dict[str, float]] = None) -> "PositionSpecificScoringMatrix":
        """Calculate log-odds scores to create a PSSM.
        
        Converts frequency probabilities to log-odds scores by comparing
        each frequency to the background frequency for that character.
        
        Args:
            background: Optional dictionary of background frequencies for
                       each alphabet character. If None, assumes uniform
                       background (equal probability for all characters).
        
        Returns:
            PositionSpecificScoringMatrix: Matrix containing log2-odds scores
                                          for sequence evaluation.
        
        Note:
            Log-odds score = log2(frequency / background_frequency)
            Positive scores indicate above-background frequency.
            Negative scores indicate below-background frequency.
        """
        if background is None:
            # Use uniform background
            background = {letter: 1.0 / len(self.alphabet) for letter in self.alphabet}
        
        pssm_values = {}
        for letter in self.alphabet:
            pssm_values[letter] = []
            for i in range(self.length):
                frequency = max(self[letter][i], MIN_PROBABILITY)  # Avoid log(0)
                bg_freq = max(background[letter], MIN_PROBABILITY)  # Avoid division by 0
                
                if DEFAULT_LOG_BASE == 2:
                    log_odds = math.log2(frequency / bg_freq)
                else:
                    log_odds = math.log(frequency / bg_freq) / math.log(DEFAULT_LOG_BASE)
                
                pssm_values[letter].append(log_odds)
        
        return PositionSpecificScoringMatrix(self.alphabet, pssm_values)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """A position-specific scoring matrix (PSSM).
    
    Contains log-odds scores for each alphabet character at each position.
    Used for scoring sequences and finding motif matches. Higher scores
    indicate better matches to the motif pattern.
    
    PSSMs are the final form of motif representation used for practical
    sequence analysis and motif searching applications.
    """

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the PositionSpecificScoringMatrix.
        
        Args:
            alphabet: String of valid characters (e.g., "ACGT").
            values: Dictionary mapping alphabet characters to lists of
                   log-odds scores. Values can be negative, positive, or zero.
        """
        super().__init__(alphabet, values)

    @property
    def min(self) -> float:
        """Return the minimum possible score.
        
        Calculates the lowest score achievable by selecting the worst
        (most negative) character at each position.
        
        Returns:
            float: Minimum possible total score for any sequence.
        """
        min_score = 0.0
        for i in range(self.length):
            position_min = min(self[letter][i] for letter in self.alphabet)
            min_score += position_min
        return min_score

    @property
    def max(self) -> float:
        """Return the maximum possible score.
        
        Calculates the highest score achievable by selecting the best
        (most positive) character at each position.
        
        Returns:
            float: Maximum possible total score for any sequence.
        """
        max_score = 0.0
        for i in range(self.length):
            position_max = max(self[letter][i] for letter in self.alphabet)
            max_score += position_max
        return max_score

    def calculate(self, sequence: str, start: int = 0) -> float:
        """Calculate the score for a sequence at the given start position.
        
        Computes the sum of log-odds scores for each position in the motif
        when aligned with the sequence starting at the specified position.
        
        Args:
            sequence: Target sequence string to score.
            start: Starting position in the sequence (0-based).
        
        Returns:
            float: Total log-odds score for the motif match.
                  Higher scores indicate better matches.
        
        Raises:
            ValueError: If the sequence is too short for the motif
                       starting at the given position.
            TypeError: If sequence is not a string or start is not an integer.
        
        Example:
            >>> pssm.calculate("ATCGATCG", start=2)
            4.52
        """
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        
        if not isinstance(start, int):
            raise TypeError("start must be an integer")
        
        if start < 0:
            raise ValueError("start position cannot be negative")
        
        if not sequence:
            raise ValueError("sequence cannot be empty")
        
        if start + self.length > len(sequence):
            raise ValueError(f"Sequence too short for motif: need {self.length} characters "
                           f"starting at position {start}, but sequence has only {len(sequence)} characters")
        
        score = 0.0
        for i in range(self.length):
            nucleotide = sequence[start + i].upper()
            if nucleotide in self.alphabet:
                score += self[nucleotide][i]
            else:
                # Handle ambiguous nucleotides by using minimum score
                min_score = min(self[letter][i] for letter in self.alphabet)
                score += min_score
        
        return score

    def search(self, sequence: str, threshold: float = 0.0) -> List[tuple]:
        """Search for motif matches in a sequence above the threshold.
        
        Slides the motif across the entire sequence and returns all
        positions where the score exceeds the specified threshold.
        
        Args:
            sequence: Target sequence string to search.
            threshold: Minimum score required for a match (default: 0.0).
        
        Returns:
            List of tuples (start, end, score) for each match above threshold.
            Positions are 0-based. 'end' is exclusive (Python slice style).
        
        Raises:
            TypeError: If sequence is not a string or threshold is not numeric.
            ValueError: If sequence is too short for the motif.
        
        Example:
            >>> matches = pssm.search("ATCGATCG", threshold=2.0)
            >>> print(matches)
            [(1, 5, 3.24), (4, 8, 2.15)]
        """
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        
        if not isinstance(threshold, (int, float)):
            raise TypeError("threshold must be a number")
        
        if not sequence:
            return []  # Empty sequence, no matches
        
        seq_len = len(sequence)
        if seq_len < self.length:
            raise ValueError(f"Sequence too short: need at least {self.length} characters, "
                           f"got {seq_len}")
        
        matches = []
        for start in range(seq_len - self.length + 1):
            try:
                score = self.calculate(sequence, start)
                if score >= threshold:
                    end = start + self.length
                    matches.append((start, end, score))
            except Exception as e:
                # Log the error but continue searching
                # In a real implementation, you might want to use proper logging
                continue
        
        return matches

    def reverse_complement(self) -> "PositionSpecificScoringMatrix":
        """Return the reverse complement of the PSSM.
        
        Creates a new PSSM representing the reverse complement motif.
        Only meaningful for DNA alphabets.
        
        Returns:
            PositionSpecificScoringMatrix: New PSSM with reverse complement
                                          pattern. Positions are reversed and
                                          complementary characters are swapped.
        
        Raises:
            ValueError: If the alphabet doesn't contain standard DNA characters.
        
        Note:
            This assumes DNA alphabet with complement mapping:
            A ↔ T, C ↔ G
        """
        # Validate that this is a DNA alphabet
        dna_chars = set(DNA_COMPLEMENT_MAP.keys())
        alphabet_chars = set(self.alphabet)
        
        if not dna_chars.issubset(alphabet_chars):
            raise ValueError(f"Reverse complement requires DNA alphabet. "
                           f"Missing characters: {dna_chars - alphabet_chars}")
        
        rc_values = {}
        for letter in self.alphabet:
            if letter in DNA_COMPLEMENT_MAP:
                complement_letter = DNA_COMPLEMENT_MAP[letter]
                rc_values[complement_letter] = list(reversed(self[letter]))
            else:
                # For non-standard letters, just reverse (maintain original behavior)
                rc_values[letter] = list(reversed(self[letter]))
        
        return PositionSpecificScoringMatrix(self.alphabet, rc_values)

    def mean(self, background: Optional[Dict[str, float]] = None) -> float:
        """Calculate the mean score of the PSSM.
        
        Computes the expected score when the PSSM is applied to random
        sequences with the given background composition.
        
        Args:
            background: Background frequency distribution. If None,
                       assumes uniform distribution.
        
        Returns:
            float: Expected score for random sequences.
                  Should be close to 0 for well-designed PSSMs.
        
        Note:
            Mean score = sum over positions of sum over characters of
            (background_freq * score)
        """
        if background is None:
            # Use uniform background
            uniform_freq = 1.0 / len(self.alphabet)
            background = {letter: uniform_freq for letter in self.alphabet}
        
        # Validate background frequencies sum to 1
        total_bg = sum(background.values())
        if abs(total_bg - 1.0) > MIN_PROBABILITY:
            # Normalize background if needed
            background = {letter: freq / total_bg for letter, freq in background.items()}
        
        mean_score = 0.0
        for i in range(self.length):
            position_mean = sum(background.get(letter, 0.0) * self[letter][i] 
                              for letter in self.alphabet)
            mean_score += position_mean
        
        return mean_score