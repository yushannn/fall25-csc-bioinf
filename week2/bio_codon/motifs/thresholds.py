"""Approximate calculation of appropriate thresholds for motif finding (Codon port)."""

from typing import Dict, List, Optional, Union
import math


class ScoreDistribution:
    """Class representing approximate score distribution for a given motif.

    Utilizes a dynamic programming approach to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurrences.
    """

    def __init__(self, motif=None, precision: int = 1000, pssm=None, background: Optional[Dict[str, float]] = None):
        """Initialize the class."""
        if pssm is None and motif is None:
            raise ValueError("Either motif or pssm must be provided")
        
        if pssm is None:
            # Use motif to get PSSM
            pssm = motif.pssm
            self.min_score = pssm.min
            self.max_score = pssm.max
            self.length = motif.length
            if background is None:
                background = motif.background
        else:
            self.min_score = pssm.min
            self.max_score = pssm.max  
            self.length = pssm.length
        
        self.interval = self.max_score - self.min_score
        self.n_points = precision * self.length
        self.step = self.interval / (self.n_points - 1) if self.n_points > 1 else 1.0
        
        # Initialize probability distributions
        self.mo_density = [0.0] * self.n_points
        self.bg_density = [0.0] * self.n_points
        
        # Set initial probability at minimum score
        min_index = self._score_to_index(self.min_score)
        if 0 <= min_index < self.n_points:
            self.mo_density[min_index] = 1.0
            self.bg_density[min_index] = 1.0
        
        # Build score distribution using dynamic programming
        self._build_distribution(pssm, background)

    def _score_to_index(self, score: float) -> int:
        """Convert score to array index."""
        if self.step == 0:
            return 0
        index = int((score - self.min_score) / self.step)
        return max(0, min(index, self.n_points - 1))

    def _index_to_score(self, index: int) -> float:
        """Convert array index to score."""
        return self.min_score + index * self.step

    def _build_distribution(self, pssm, background: Dict[str, float]) -> None:
        """Build score distribution using dynamic programming."""
        if background is None:
            # Use uniform background
            background = {letter: 1.0 / len(pssm.alphabet) for letter in pssm.alphabet}
        
        # For each position in the motif
        for position in range(self.length):
            mo_new = [0.0] * self.n_points
            bg_new = [0.0] * self.n_points
            
            # For each possible nucleotide at this position
            for letter in pssm.alphabet:
                score_contribution = pssm[letter][position]
                bg_prob = background[letter]
                
                # Calculate motif probability (assuming position weight matrix exists)
                mo_prob = bg_prob * (2 ** score_contribution)  # Approximate
                
                # Update distributions
                score_index_diff = self._score_to_index(score_contribution) - self._score_to_index(0)
                
                for i in range(self.n_points):
                    new_index = i + score_index_diff
                    if 0 <= new_index < self.n_points:
                        mo_new[new_index] += self.mo_density[i] * mo_prob
                        bg_new[new_index] += self.bg_density[i] * bg_prob
            
            self.mo_density = mo_new
            self.bg_density = bg_new

    def threshold_fpr(self, fpr: float) -> float:
        """Calculate threshold for a given false positive rate."""
        if not (0 <= fpr <= 1):
            raise ValueError("FPR must be between 0 and 1")
        
        # Calculate cumulative background distribution (from high to low scores)
        cumulative_bg = 0.0
        total_bg = sum(self.bg_density)
        
        if total_bg == 0:
            return self.max_score
        
        for i in range(self.n_points - 1, -1, -1):
            cumulative_bg += self.bg_density[i]
            current_fpr = cumulative_bg / total_bg
            
            if current_fpr >= fpr:
                return self._index_to_score(i)
        
        return self.min_score

    def threshold_fnr(self, fnr: float) -> float:
        """Calculate threshold for a given false negative rate."""
        if not (0 <= fnr <= 1):
            raise ValueError("FNR must be between 0 and 1")
        
        # Calculate cumulative motif distribution (from low to high scores)
        cumulative_mo = 0.0
        total_mo = sum(self.mo_density)
        
        if total_mo == 0:
            return self.min_score
        
        for i in range(self.n_points):
            cumulative_mo += self.mo_density[i]
            current_fnr = cumulative_mo / total_mo
            
            if current_fnr >= fnr:
                return self._index_to_score(i)
        
        return self.max_score

    def threshold_balanced(self) -> float:
        """Calculate threshold that balances false positive and false negative rates."""
        best_threshold = self.min_score
        best_difference = float('inf')
        
        total_mo = sum(self.mo_density)
        total_bg = sum(self.bg_density)
        
        if total_mo == 0 or total_bg == 0:
            return (self.min_score + self.max_score) / 2
        
        cumulative_mo = 0.0
        cumulative_bg = sum(self.bg_density)  # Start from total (high scores)
        
        for i in range(self.n_points):
            # Update cumulative distributions
            cumulative_mo += self.mo_density[i]
            if i > 0:
                cumulative_bg -= self.bg_density[self.n_points - i]
            
            # Calculate rates
            fnr = cumulative_mo / total_mo
            fpr = cumulative_bg / total_bg
            
            # Find point where FPR ≈ FNR
            difference = abs(fpr - fnr)
            if difference < best_difference:
                best_difference = difference
                best_threshold = self._index_to_score(i)
        
        return best_threshold

    def threshold_patser(self) -> float:
        """Calculate PATSER-style threshold (information content based)."""
        # Simple implementation: use 2/3 of maximum possible score
        return self.min_score + 0.67 * (self.max_score - self.min_score)

    def score_to_fpr(self, score: float) -> float:
        """Calculate false positive rate for a given score threshold."""
        total_bg = sum(self.bg_density)
        if total_bg == 0:
            return 0.0
        
        threshold_index = self._score_to_index(score)
        cumulative_bg = sum(self.bg_density[threshold_index:])
        
        return cumulative_bg / total_bg

    def score_to_fnr(self, score: float) -> float:
        """Calculate false negative rate for a given score threshold."""
        total_mo = sum(self.mo_density)
        if total_mo == 0:
            return 1.0
        
        threshold_index = self._score_to_index(score)
        cumulative_mo = sum(self.mo_density[:threshold_index])
        
        return cumulative_mo / total_mo


def score_distribution(motif, background: Optional[Dict[str, float]] = None, precision: int = 1000) -> ScoreDistribution:
    """Create a ScoreDistribution object for the given motif."""
    return ScoreDistribution(motif=motif, background=background, precision=precision)


def threshold_from_pvalue(pssm, pvalue: float, background: Optional[Dict[str, float]] = None, precision: int = 1000) -> float:
    """Calculate threshold from p-value using score distribution."""
    distribution = ScoreDistribution(pssm=pssm, background=background, precision=precision)
    return distribution.threshold_fpr(pvalue)


def pvalue_from_threshold(pssm, threshold: float, background: Optional[Dict[str, float]] = None, precision: int = 1000) -> float:
    """Calculate p-value from threshold using score distribution."""
    distribution = ScoreDistribution(pssm=pssm, background=background, precision=precision)
    return distribution.score_to_fpr(threshold)