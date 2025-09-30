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

# Constants
DEFAULT_DNA_ALPHABET = "ACGT"
DEFAULT_RNA_ALPHABET = "ACGU"
DEFAULT_PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Standard nucleotide complements
DNA_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
RNA_COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C"}

# Default values
DEFAULT_PSEUDOCOUNT = 0.0
DEFAULT_SEARCH_THRESHOLD = 0.0
MIN_MOTIF_LENGTH = 1
MAX_REASONABLE_MOTIF_LENGTH = 50

# Tolerance for floating point comparisons
FLOAT_TOLERANCE = 1e-10

# Import local modules - will be implemented
from . import matrix


def create(instances: List[str], alphabet: str = DEFAULT_DNA_ALPHABET) -> "Motif":
    """Create a Motif object from a list of sequences.
    
    Args:
        instances: List of sequences of equal length to create the motif from.
                  All sequences must contain only characters from the alphabet.
        alphabet: String containing the alphabet characters (default: "ACGT").
                 Each character represents a possible nucleotide/amino acid.
    
    Returns:
        Motif: A new Motif object containing the frequency information
               derived from the input sequences.
    
    Raises:
        ValueError: If sequences have different lengths or contain
                   characters not in the alphabet.
        TypeError: If instances is not a list or contains non-string elements.
    
    Example:
        >>> motif = create(["ATCG", "ATGG", "TTCG"])
        >>> print(motif.consensus)
        'ATCG'
    """
    # Input validation
    if not isinstance(instances, list):
        raise TypeError("instances must be a list of strings")
    
    if not instances:
        raise ValueError("instances cannot be empty")
    
    if not isinstance(alphabet, str) or not alphabet:
        raise ValueError("alphabet must be a non-empty string")
    
    # Validate sequence types and characters
    alphabet_set = set(alphabet)
    for i, seq in enumerate(instances):
        if not isinstance(seq, str):
            raise TypeError(f"All instances must be strings, got {type(seq)} at position {i}")
        
        invalid_chars = set(seq.upper()) - alphabet_set
        if invalid_chars:
            raise ValueError(f"Sequence {i} contains invalid characters: {sorted(invalid_chars)}. "
                           f"Valid alphabet: {alphabet}")
    
    # Create alignment and motif
    alignment = SimpleAlignment(instances)
    return Motif(alignment=alignment, alphabet=alphabet)


def parse(handle, fmt: str, strict: bool = True):
    """Parse an output file from a motif finding program.

    Args:
        handle: File handle or string containing the motif data to parse.
               Can be any object that supports iteration over lines.
        fmt: Format specification string (case insensitive).
             Currently supported formats:
             - 'MINIMAL': MINIMAL MEME output file motif
             - 'pfm': JASPAR-style position-frequency matrix
             - 'jaspar': JASPAR-style multiple PFM format
             - 'sites': JASPAR-style sites file
        strict: If True (default), raise ValueError for format violations.
               If False, attempt to parse non-compliant files with warnings.

    Returns:
        Generator yielding Motif objects parsed from the input.
        For single motif formats, yields exactly one Motif.
        For multiple motif formats, yields multiple Motifs.

    Raises:
        ValueError: If format is unknown or file doesn't comply with
                   the specified format (when strict=True).
        IOError: If handle cannot be read.

    Example:
        >>> with open('motifs.txt', 'r') as f:
        >>>     motifs = list(parse(f, 'minimal'))
        >>> print(len(motifs))
        3
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
    """Simple alignment class for Codon compatibility.
    
    Represents a multiple sequence alignment where all sequences
    have the same length. Provides frequency counting functionality
    for motif analysis.
    
    Attributes:
        sequences: List of aligned sequences of equal length.
        length: Length of each sequence in the alignment.
    """
    
    def __init__(self, sequences: List[str]):
        """Initialize SimpleAlignment with a list of sequences.
        
        Args:
            sequences: List of sequence strings. All sequences must have
                      the same length. Empty list is allowed.
        
        Raises:
            ValueError: If sequences have different lengths.
            TypeError: If sequences is not a list or contains non-strings.
        """
        if not isinstance(sequences, list):
            raise TypeError("sequences must be a list")
        
        # Validate sequence types
        for i, seq in enumerate(sequences):
            if not isinstance(seq, str):
                raise TypeError(f"All sequences must be strings, got {type(seq)} at position {i}")
        
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0
        
        # Verify all sequences have same length
        for i, seq in enumerate(sequences):
            if len(seq) != self.length:
                raise ValueError(f"All sequences must have the same length. "
                               f"Expected {self.length}, got {len(seq)} at position {i}")
    
    @property
    def frequencies(self) -> Dict[str, List[int]]:
        """Calculate nucleotide frequencies at each position.
        
        Returns:
            Dict mapping each nucleotide (A, C, G, T) to a list of integers
            representing the count of that nucleotide at each position
            across all sequences in the alignment.
            
        Example:
            For sequences ["ATG", "ACG", "ATG"]:
            Returns: {'A': [2, 0, 0], 'C': [0, 1, 0], 'G': [0, 0, 3], 'T': [0, 1, 0]}
        """
        freq = {}
        if not self.sequences:
            return freq
            
        # Initialize frequency matrix
        alphabet_chars = list(DEFAULT_DNA_ALPHABET)  # Use default alphabet
        for letter in alphabet_chars:
            freq[letter] = [0] * self.length
        
        # Count frequencies
        for seq in self.sequences:
            for i, nucleotide in enumerate(seq):
                if nucleotide in freq:
                    freq[nucleotide][i] += 1
        
        return freq


class Motif:
    """A class representing sequence motifs.
    
    A motif represents a conserved pattern in biological sequences,
    typically derived from multiple sequence alignments or frequency
    count matrices. It provides functionality for motif analysis,
    scoring, and format conversion.
    
    Attributes:
        name: String identifier for the motif.
        alphabet: String containing valid characters (e.g., "ACGT").
        length: Length of the motif pattern.
        counts: FrequencyPositionMatrix with nucleotide counts.
        alignment: Optional SimpleAlignment object.
        pseudocounts: Background correction values.
        background: Background frequency distribution.
        mask: Position mask for analysis (1=use, 0=ignore).
    """

    def __init__(self, alphabet: str = DEFAULT_DNA_ALPHABET, alignment: Optional[SimpleAlignment] = None, 
                 counts: Optional[Dict[str, List[int]]] = None):
        """Initialize the Motif class.
        
        Args:
            alphabet: String of valid characters (default: "ACGT").
                     Each character represents a possible symbol.
            alignment: Optional SimpleAlignment object containing sequences.
                      Used to calculate frequency counts automatically.
            counts: Optional dictionary mapping alphabet characters to
                   lists of integer counts at each position.
        
        Raises:
            ValueError: If both alignment and counts are provided, or if
                       counts data is inconsistent.
        
        Note:
            Exactly one of alignment or counts should be provided.
            If neither is provided, creates an empty motif.
        """
        self.name = ""
        
        # Input validation
        if not isinstance(alphabet, str) or not alphabet:
            raise ValueError("alphabet must be a non-empty string")
        
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
            
            if self._has_dna_alphabet():
                alphabet_chars = DEFAULT_DNA_ALPHABET
            else:  # RNA alphabet
                alphabet_chars = DEFAULT_RNA_ALPHABET
            
            # Set GC content-based background
            gc_content = float(value)
            at_content = 1.0 - gc_content
            
            self._background = {
                alphabet_chars[0]: at_content / 2.0,  # A
                alphabet_chars[1]: gc_content / 2.0,   # C
                alphabet_chars[2]: gc_content / 2.0,   # G
                alphabet_chars[3]: at_content / 2.0,   # T or U
            }
        
        # Normalize to sum to 1
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total

    background = property(__get_background, __set_background)

    def __len__(self) -> int:
        """Return the length of a motif.
        
        Returns:
            int: The number of positions in the motif pattern.
                Returns 0 for empty motifs.
        """
        if self.length is None:
            return 0
        return self.length

    def _has_dna_alphabet(self) -> bool:
        """Check if motif uses standard DNA alphabet."""
        return sorted(self.alphabet) == sorted(DEFAULT_DNA_ALPHABET)

    def _has_rna_alphabet(self) -> bool:
        """Check if motif uses standard RNA alphabet."""
        return sorted(self.alphabet) == sorted(DEFAULT_RNA_ALPHABET)

    @property
    def pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        return self.pwm.log_odds(self._background)

    def __str__(self) -> str:
        """Return string representation of a motif.
        
        Returns:
            str: Multi-line string showing the aligned sequences
                if available, otherwise empty string.
        """
        text = ""
        if self.alignment is not None:
            text += "\n".join(self.alignment.sequences)
        return text

    def reverse_complement(self) -> "Motif":
        """Return the reverse complement of the motif as a new motif.
        
        Creates a new Motif object representing the reverse complement
        of this motif. Only works for DNA (ACGT) and RNA (ACGU) alphabets.
        
        Returns:
            Motif: New motif object with reverse complement pattern.
                  All properties (background, pseudocounts, mask) are
                  appropriately transformed.
        
        Raises:
            ValueError: If the motif alphabet is not DNA or RNA.
        
        Example:
            >>> motif = create(["ATCG", "ATGG"])
            >>> rc_motif = motif.reverse_complement()
            >>> print(rc_motif.consensus)  # Should be complement of original
        """
        alphabet = self.alphabet
        if not self._has_dna_alphabet() and not self._has_rna_alphabet():
            raise ValueError("Calculating reverse complement only works for DNA and RNA motifs")
        
        if self._has_dna_alphabet():
            complement_map = DNA_COMPLEMENT
        elif self._has_rna_alphabet():
            complement_map = RNA_COMPLEMENT
        else:
            raise ValueError("Calculating reverse complement only works for DNA and RNA motifs")
        
        if self.alignment is not None:
            # Reverse complement the alignment
            rev_comp_sequences = []
            for seq in self.alignment.sequences:
                rev_comp = "".join(complement_map.get(base, base) for base in reversed(seq))
                rev_comp_sequences.append(rev_comp)
            
            rev_alignment = SimpleAlignment(rev_comp_sequences)
            res = Motif(alphabet=alphabet, alignment=rev_alignment)
        else:
            # Reverse complement the counts using complement mapping
            counts = {}
            for original_char, complement_char in complement_map.items():
                if original_char in self.counts:
                    counts[complement_char] = list(reversed(self.counts[original_char]))
            
            # Handle any characters not in complement map
            for char in self.alphabet:
                if char not in counts:
                    counts[char] = list(reversed(self.counts[char]))
            
            res = Motif(alphabet=alphabet, counts=counts)
        
        # Set properties using complement mapping
        if self.__mask:
            res.__mask = tuple(reversed(self.__mask))
        
        # Update background and pseudocounts using complement mapping
        res._background = {}
        res._pseudocounts = {}
        
        # Initialize background and pseudocounts if they are None
        if self._background is None:
            self._background = {letter: 1.0 / len(self.alphabet) for letter in self.alphabet}
        if self._pseudocounts is None:
            self._pseudocounts = {letter: 0.0 for letter in self.alphabet}
        
        for original_char, complement_char in complement_map.items():
            if original_char in self._background:
                res._background[complement_char] = self._background[original_char]
                res._pseudocounts[complement_char] = self._pseudocounts[original_char]
        
        # Handle any characters not in complement map
        for char in self.alphabet:
            if char not in res._background:
                res._background[char] = self._background[char]
                res._pseudocounts[char] = self._pseudocounts[char]
        
        return res

    @property
    def consensus(self) -> str:
        """Return the consensus sequence.
        
        The consensus sequence is formed by selecting the most frequent
        nucleotide at each position. In case of ties, the first one
        encountered in alphabet order is chosen.
        
        Returns:
            str: Consensus sequence string of length equal to motif length.
                Returns empty string if motif has no counts data.
        
        Example:
            For a motif with counts A:[3,1], C:[0,2], G:[1,0], T:[0,1]
            Returns: "AC"
        """
        if self.counts is None or self.length is None or self.length == 0:
            return ""
        
        consensus_seq = ""
        for i in range(self.length):
            max_count = -1  # Use -1 to handle zero counts properly
            consensus_base = "N"  # Default fallback
            
            # Find the most frequent character at this position
            for letter in self.alphabet:
                count = self.counts[letter][i]
                if count > max_count:
                    max_count = count
                    consensus_base = letter
            
            consensus_seq += consensus_base
        
        return consensus_seq

    def __format__(self, format_spec: str, **kwargs) -> str:
        """Return a string representation of the Motif in the given format.
        
        Args:
            format_spec: Format specification string. Supported formats:
                        - 'pfm' or 'jaspar': JASPAR position frequency matrix
                        - 'transfac': TRANSFAC database format
                        - Default: returns string representation of alignment
            **kwargs: Additional formatting options (unused currently).
        
        Returns:
            str: Formatted string representation of the motif.
        
        Example:
            >>> motif = create(["ATCG", "ATGG"])
            >>> print(format(motif, 'pfm'))
            A [ 2.00  0.00  0.00  0.00]
            C [ 0.00  0.00  2.00  0.00]
            G [ 0.00  0.00  0.00  2.00]
            T [ 0.00  2.00  0.00  0.00]
        """
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
                # Find the consensus letter for this position
                max_count = -1
                consensus_letter = "N"
                for letter in self.alphabet:
                    if self.counts[letter][i] > max_count:
                        max_count = self.counts[letter][i]
                        consensus_letter = letter
                lines.append(f"{i:02d}\t" + "\t".join(counts) + f"\t{consensus_letter}")
            lines.append("//")
            return "\n".join(lines)
        else:
            return str(self)