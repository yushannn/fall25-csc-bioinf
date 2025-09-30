"""Tools for sequence motif analysis (Codon port).

Bio.motifs contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also includes functionality for parsing output from the AlignACE, MEME,
and MAST programs, as well as files in the TRANSFAC format.
"""

# Import required modules
# Python imports will be handled conditionally in test file
# For now, implement core functionality without external dependencies
from typing import Dict, List, Optional, Union, Any, Tuple

# Import local modules - will be implemented
from . import matrix


def create(instances: List[str], alphabet: str = "ACGT") -> "Motif":
    """Create a Motif object from a list of sequences."""
    # For now, create a simple alignment-like structure
    alignment = SimpleAlignment(instances)
    return Motif(alignment=alignment, alphabet=alphabet)


def parse(handle, fmt: str, strict: bool = True):
    """Parse an output file from a motif finding program.

    Currently supported formats (case is ignored):
     - MINIMAL:          MINIMAL MEME output file motif
     - pfm:              JASPAR-style position-frequency matrix
     - jaspar:           JASPAR-style multiple PFM format
     - sites:            JASPAR-style sites file

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    
    if fmt == "minimal":
        from . import minimal
        return minimal.read(handle)
    # Add other formats as needed
    else:
        raise ValueError(f"Unknown format: {fmt}")


def read(handle, fmt: str, strict: bool = True):
    """Read a single motif from a file."""
    motifs = list(parse(handle, fmt, strict))
    if len(motifs) == 0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    motif = motifs[0]
    return motif


class SimpleAlignment:
    """Simple alignment class for Codon compatibility."""
    
    def __init__(self, sequences: List[str]):
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0
        # Verify all sequences have same length
        for seq in sequences:
            if len(seq) != self.length:
                raise ValueError("All sequences must have the same length")
    
    @property
    def frequencies(self) -> Dict[str, List[int]]:
        """Calculate nucleotide frequencies at each position."""
        freq = {}
        if not self.sequences:
            return freq
            
        # Initialize frequency matrix
        alphabet = ["A", "C", "G", "T"]
        for letter in alphabet:
            freq[letter] = [0] * self.length
        
        # Count frequencies
        for seq in self.sequences:
            for i, nucleotide in enumerate(seq):
                if nucleotide in freq:
                    freq[nucleotide][i] += 1
        
        return freq


class Motif:
    """A class representing sequence motifs."""

    def __init__(self, alphabet: str = "ACGT", alignment: Optional[SimpleAlignment] = None, 
                 counts: Optional[Dict[str, List[int]]] = None):
        """Initialize the Motif class."""
        self.name = ""
        
        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            self.alignment = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet, counts)
            self.length = self.counts.length
        elif alignment is not None:
            length = alignment.length
            frequencies = alignment.frequencies
            # Ensure all alphabet letters are present
            for letter in alphabet:
                if letter not in frequencies:
                    frequencies[letter] = [0] * length
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
        else:
            self.counts = None
            self.alignment = None
            self.length = None
            
        self.alphabet = alphabet
        self._pseudocounts = None
        self._background = None
        self.__mask = None

    def __get_mask(self) -> Tuple[int, ...]:
        return self.__mask if self.__mask is not None else ()

    def __set_mask(self, mask):
        if self.length is None:
            self.__mask = ()
        elif mask is None:
            self.__mask = tuple(1 for _ in range(self.length))
        elif len(mask) != self.length:
            raise ValueError(
                f"The length ({len(mask)}) of the mask is inconsistent with the length ({self.length}) of the motif"
            )
        elif isinstance(mask, str):
            self.__mask = []
            for char in mask:
                if char == "*":
                    self.__mask.append(1)
                elif char == " ":
                    self.__mask.append(0)
                else:
                    raise ValueError(f"Mask should contain only '*' or ' ' and not a '{char}'")
            self.__mask = tuple(self.__mask)
        else:
            self.__mask = tuple(int(bool(c)) for c in mask)

    mask = property(__get_mask, __set_mask)

    def __get_pseudocounts(self) -> Dict[str, float]:
        if self._pseudocounts is None:
            return {letter: 0.0 for letter in self.alphabet}
        return self._pseudocounts

    def __set_pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, dict):
            self._pseudocounts = {letter: float(value[letter]) for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = {letter: float(value) for letter in self.alphabet}

    pseudocounts = property(__get_pseudocounts, __set_pseudocounts)

    def __get_background(self) -> Dict[str, float]:
        if self._background is None:
            return {letter: 1.0 for letter in self.alphabet}
        return self._background

    def __set_background(self, value):
        if isinstance(value, dict):
            self._background = {letter: float(value[letter]) for letter in self.alphabet}
        elif value is None:
            self._background = {letter: 1.0 for letter in self.alphabet}
        else:
            if not self._has_dna_alphabet() and not self._has_rna_alphabet():
                raise ValueError(
                    "Setting the background to a single value only works for DNA and RNA motifs"
                )
            T_or_U = "T" if self._has_dna_alphabet() else "U"
            self._background = {}
            self._background["A"] = (1.0 - value) / 2.0
            self._background["C"] = value / 2.0
            self._background["G"] = value / 2.0
            self._background[T_or_U] = (1.0 - value) / 2.0
        
        # Normalize to sum to 1
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total

    background = property(__get_background, __set_background)

    def __len__(self) -> int:
        """Return the length of a motif."""
        if self.length is None:
            return 0
        return self.length

    def _has_dna_alphabet(self) -> bool:
        return sorted(self.alphabet) == ["A", "C", "G", "T"]

    def _has_rna_alphabet(self) -> bool:
        return sorted(self.alphabet) == ["A", "C", "G", "U"]

    @property
    def pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        return self.pwm.log_odds(self._background)

    def __str__(self) -> str:
        """Return string representation of a motif."""
        text = ""
        if self.alignment is not None:
            text += "\n".join(self.alignment.sequences)
        return text

    def reverse_complement(self) -> "Motif":
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if not self._has_dna_alphabet() and not self._has_rna_alphabet():
            raise ValueError("Calculating reverse complement only works for DNA and RNA motifs")
        
        T_or_U = "T" if self._has_dna_alphabet() else "U"
        
        if self.alignment is not None:
            # Reverse complement the alignment
            rev_comp_sequences = []
            complement_map = {"A": T_or_U, "C": "G", "G": "C", T_or_U: "A"}
            for seq in self.alignment.sequences:
                rev_comp = "".join(complement_map.get(base, base) for base in reversed(seq))
                rev_comp_sequences.append(rev_comp)
            
            rev_alignment = SimpleAlignment(rev_comp_sequences)
            res = Motif(alphabet=alphabet, alignment=rev_alignment)
        else:
            # Reverse complement the counts
            counts = {
                "A": list(reversed(self.counts[T_or_U])),
                "C": list(reversed(self.counts["G"])),
                "G": list(reversed(self.counts["C"])),
                T_or_U: list(reversed(self.counts["A"])),
            }
            res = Motif(alphabet=alphabet, counts=counts)
        
        # Set properties
        if self.__mask:
            res.__mask = tuple(reversed(self.__mask))
        res._background = {
            "A": self._background[T_or_U],
            "C": self._background["G"],
            "G": self._background["C"],
            T_or_U: self._background["A"],
        }
        res._pseudocounts = {
            "A": self._pseudocounts[T_or_U],
            "C": self._pseudocounts["G"],
            "G": self._pseudocounts["C"],
            T_or_U: self._pseudocounts["A"],
        }
        
        return res

    @property
    def consensus(self) -> str:
        """Return the consensus sequence."""
        if self.counts is None:
            return ""
        
        consensus_seq = ""
        for i in range(self.length):
            max_count = 0
            consensus_base = "N"
            for letter in self.alphabet:
                if self.counts[letter][i] > max_count:
                    max_count = self.counts[letter][i]
                    consensus_base = letter
            consensus_seq += consensus_base
        
        return consensus_seq

    def __format__(self, format_spec: str, **kwargs) -> str:
        """Return a string representation of the Motif in the given format."""
        if format_spec in ("pfm", "jaspar"):
            # Simple PFM format implementation
            lines = []
            for letter in self.alphabet:
                counts_str = " ".join(str(count) for count in self.counts[letter])
                lines.append(f"{letter} [{counts_str}]")
            return "\n".join(lines)
        elif format_spec == "transfac":
            # Simple TRANSFAC-like format
            lines = [f"DE\t{self.name}"]
            for i in range(self.length):
                counts = [str(self.counts[letter][i]) for letter in self.alphabet]
                lines.append(f"{i:02d}\t" + "\t".join(counts) + f"\t{letter}")
            lines.append("//")
            return "\n".join(lines)
        else:
            return str(self)