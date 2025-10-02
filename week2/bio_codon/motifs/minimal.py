"""Module for the support of MEME minimal motif format (Codon port)."""

from typing import List, Dict, Optional, TextIO, Any
from . import Motif

# Constants for MEME format parsing
MEME_VERSION_PREFIX = "MEME version"
ALPHABET_PREFIX = "ALPHABET="
BACKGROUND_PREFIX = "Background letter frequencies"
MOTIF_PREFIX = "letter-probability matrix:"
MOTIF_NAME_KEYWORDS = ["MOTIF", "motif"]
DEFAULT_MOTIF_NAME = "Unnamed"
MAX_LINE_LENGTH = 10000  # Safety limit for line parsing
DEFAULT_BACKGROUND_FREQ = 0.25  # Default frequency if not specified


class Record:
    """Record holding MEME minimal format data.
    
    A Record object holds the information from a MEME minimal format file,
    including version information, alphabet, background frequencies, and
    a list of motifs.
    
    Attributes:
        version: MEME version string.
        alphabet: String containing alphabet characters.
        background: Dictionary of background character frequencies.
        motifs: List of Motif objects parsed from the file.
    """

    def __init__(self):
        """Initialize an empty MEME Record.
        
        Creates a new Record with empty fields that will be populated
        during parsing of a MEME minimal format file.
        """
        self.version: str = ""
        self.alphabet: str = ""
        self.background: Dict[str, float] = {}
        self.motifs: List[Motif] = []

    def __len__(self) -> int:
        """Return number of motifs in record."""
        return len(self.motifs)

    def __getitem__(self, key) -> Motif:
        """Get motif by index or name.
        
        Args:
            key: Either an integer index or string name.
        
        Returns:
            Motif: The requested motif.
        
        Raises:
            IndexError: If integer index is out of range.
            KeyError: If string name is not found.
            TypeError: If key is neither int nor str.
        """
        if isinstance(key, int):
            if not self.motifs:
                raise IndexError("No motifs in record")
            return self.motifs[key]  # Let list handle IndexError
        elif isinstance(key, str):
            if not key.strip():
                raise KeyError("Motif name cannot be empty")
            # Search by name
            for motif in self.motifs:
                if hasattr(motif, 'name') and motif.name == key:
                    return motif
            raise KeyError(f"Motif with name '{key}' not found")
        else:
            raise TypeError(f"Key must be int or str, got {type(key)}")

    def __iter__(self):
        """Iterate over motifs."""
        return iter(self.motifs)


def read(handle: TextIO) -> Record:
    """Parse the text output of the MEME program into a Record object.
    
    Parses MEME minimal format files containing motif data. The function
    extracts version information, alphabet, background frequencies, and
    position weight matrices for each motif.
    
    Args:
        handle: File handle or text stream to read from. Must support
               line iteration.
    
    Returns:
        Record: A Record object containing all parsed information.
    
    Raises:
        ValueError: If the file format is invalid or corrupted.
        IOError: If there are problems reading from the handle.
    
    Examples:
        >>> from bio_codon.motifs import minimal
        >>> with open("motifs/meme.out") as f:
        ...     record = minimal.read(f)
        ...
        >>> for motif in record:
        ...     print(motif.name, getattr(motif, 'evalue', 'N/A'))
    """
    record = Record()
    
    # Read version
    _read_version(record, handle)
    
    # Read alphabet  
    _read_alphabet(record, handle)
    
    # Read background
    _read_background(record, handle)
    
    # Read motifs
    motif_number = 0
    while True:
        motif = _read_motif(record, handle, motif_number)
        if motif is None:
            break
        record.motifs.append(motif)
        motif_number += 1
    
    return record


def _read_version(record: Record, handle: TextIO) -> None:
    """Read MEME version information."""
    for line in handle:
        line = line.strip()
        if line.startswith("MEME version"):
            record.version = line.split()[-1]
            break


def _read_alphabet(record: Record, handle: TextIO) -> None:
    """Read alphabet definition."""
    for line in handle:
        line = line.strip()
        if line.startswith("ALPHABET"):
            # Extract alphabet - could be DNA, RNA, or protein
            parts = line.split()
            if len(parts) > 1:
                alphabet_part = parts[1]
                if alphabet_part in ("ACGT", "DNA"):
                    record.alphabet = "ACGT" 
                elif alphabet_part in ("ACGU", "RNA"):
                    record.alphabet = "ACGU"
                else:
                    # Custom alphabet
                    record.alphabet = alphabet_part
            break
    
    # Default to DNA if not specified
    if not record.alphabet:
        record.alphabet = "ACGT"


def _read_background(record: Record, handle: TextIO) -> None:
    """Read background frequencies."""
    for line in handle:
        line = line.strip()
        if line.startswith("Background letter frequencies"):
            # Read next line with frequencies
            freq_line = next(handle, "").strip()
            if freq_line:
                # Parse frequencies like "A 0.25 C 0.25 G 0.25 T 0.25"
                parts = freq_line.split()
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts):
                        letter = parts[i]
                        freq = float(parts[i + 1])
                        record.background[letter] = freq
            break
    
    # Default uniform background if not specified
    if not record.background:
        n_letters = len(record.alphabet)
        uniform_freq = 1.0 / n_letters
        record.background = {letter: uniform_freq for letter in record.alphabet}


def _read_motif(record: Record, handle: TextIO, motif_number: int) -> Optional[Motif]:
    """Read a single motif from the file."""
    motif_found = False
    motif_name = ""
    width = 0
    nsites = 0
    evalue = 0.0
    letter_probability_matrix = []
    
    # Look for motif header
    for line in handle:
        line = line.strip()
        if line.startswith("MOTIF"):
            motif_found = True
            parts = line.split()
            if len(parts) > 1:
                motif_name = parts[1]
            else:
                motif_name = f"motif_{motif_number + 1}"
            break
    
    if not motif_found:
        return None
    
    # Read motif parameters
    for line in handle:
        line = line.strip()
        if line.startswith("letter-probability matrix"):
            # Parse parameters from this line
            # Format: "letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 1.2e-05"
            parts = line.split()
            for part in parts:
                if part.startswith("w="):
                    width = int(part[2:])
                elif part.startswith("nsites="):
                    nsites = int(part[7:])
                elif part.startswith("E="):
                    evalue = float(part[2:])
            break
    
    # Read the probability matrix
    for i in range(width):
        line = next(handle, "").strip()
        if line:
            probs = [float(x) for x in line.split()]
            if len(probs) == len(record.alphabet):
                letter_probability_matrix.append(probs)
    
    # Convert matrix format (position x letter) to (letter x position)
    if letter_probability_matrix:
        pwm_dict = {}
        for j, letter in enumerate(record.alphabet):
            pwm_dict[letter] = [letter_probability_matrix[i][j] for i in range(width)]
        
        # Create motif with PWM data
        from . import matrix
        counts_dict = {}
        # Convert probabilities to approximate counts
        for letter in record.alphabet:
            counts_dict[letter] = [int(prob * nsites) for prob in pwm_dict[letter]]
        
        motif = Motif(alphabet=record.alphabet, counts=counts_dict)
        motif.name = motif_name
        # Store additional attributes
        motif.evalue = evalue
        motif.nsites = nsites
        
        return motif
    
    return None


def write(motifs: List[Motif]) -> str:
    """Write motifs in MEME minimal format."""
    lines = []
    
    # Header
    lines.append("MEME version 4.0 (Codon port)")
    lines.append("")
    
    # Alphabet (assume DNA for now)
    lines.append("ALPHABET= ACGT")
    lines.append("")
    
    # Background
    lines.append("Background letter frequencies")
    bg_line = ""
    for letter in "ACGT":
        bg_line += f"{letter} 0.25 "
    lines.append(bg_line.strip())
    lines.append("")
    
    # Motifs
    for i, motif in enumerate(motifs):
        lines.append(f"MOTIF {motif.name or f'motif_{i+1}'}")
        lines.append("")
        
        # Calculate sites count (approximate)
        if motif.counts:
            nsites = max(sum(motif.counts[letter]) for letter in motif.alphabet) // motif.length
        else:
            nsites = 10  # Default
        
        lines.append(f"letter-probability matrix: alength= {len(motif.alphabet)} w= {motif.length} nsites= {nsites} E= 1e-06")
        
        # Write probability matrix
        pwm = motif.pwm
        for i in range(motif.length):
            prob_line = ""
            for letter in motif.alphabet:
                prob_line += f"{pwm[letter][i]:.6f} "
            lines.append(prob_line.strip())
        
        lines.append("")
    
    return "\n".join(lines)