import os
from abc import ABC, abstractmethod
from typing import Union, Tuple

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.

    Defines the common interface for DNA, RNA, and protein sequences,
    including length support, indexing, slicing, string representation,
    and alphabet validation.
    """

    def __init__(self, sequence: str):
        """
        Initialize a biological sequence.

        Arguments:
            sequence: The raw sequence string (case-insensitive).
        """
        self._sequence = sequence.upper().strip()

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self._sequence)

    def __getitem__(self, index: Union[int, slice]) -> str:
        """
        Support indexing and slicing.

        Returns:
            A substring or single character.
        """
        return self._sequence[index]

    def __str__(self) -> str:
        """Return a human-readable string representation."""
        return self._sequence

    def __repr__(self) -> str:
        """Return a developer-friendly representation."""
        return f"{self.__class__.__name__}('{self._sequence}')"

    @abstractmethod
    def check_alphabet(self) -> bool:
        """
        Validate that the sequence contains only valid characters
        for this type of biological sequence.

        Returns:
            True if the sequence is valid, False otherwise.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for nucleic acid sequences (DNA and RNA).

    Implements complement, reverse, and reverse_complement methods
    using polymorphic class attributes for base mapping and alphabet validation.
    """

    _complement_map: dict[str, str] = {}
    _valid_bases: set[str] = set()

    def __init__(self, sequence: str):
        """
        Initialize a nucleic acid sequence.

        Arguments:
            sequence: The raw nucleotide sequence string.

        Raises:
            NotImplementedError: If instantiated directly from NucleicAcidSequence.
            ValueError: If sequence contains invalid characters.
        """
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("NucleicAcidSequence is an abstract base class")
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError(f"Invalid characters for {self.__class__.__name__}: {self._sequence}")

    def check_alphabet(self) -> bool:
        """Validate that sequence contains only valid nucleotide bases."""
        return set(self._sequence).issubset(self.__class__._valid_bases)

    def complement(self) -> "NucleicAcidSequence":
        """
        Return the complementary sequence.

        Uses polymorphic _complement_map defined in subclasses.
        Returns an instance of the same class as self.
        """
        complemented = "".join(
            self.__class__._complement_map.get(base, base) for base in self._sequence
        )
        return self.__class__(complemented)

    def reverse(self) -> "NucleicAcidSequence":
        """Return the reversed sequence."""
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """Return the reverse complement of the sequence."""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    Class representing a DNA sequence.

    Implements DNA-specific alphabet validation and complement mapping.
    Provides transcribe method to convert DNA to RNA.
    """

    _complement_map: dict[str, str] = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    _valid_bases: set[str] = {'A', 'T', 'G', 'C'}

    def transcribe(self) -> "RNASequence":
        """
        Transcribe DNA to RNA by replacing T with U.

        Returns:
            RNASequence object with transcribed sequence.
        """
        rna_seq = self._sequence.replace('T', 'U')
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    """
    Class representing an RNA sequence.

    Implements RNA-specific alphabet validation and complement mapping.
    """

    _complement_map: dict[str, str] = {
        'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
        'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'
    }
    _valid_bases: set[str] = {'A', 'U', 'G', 'C'}


class AminoAcidSequence(BiologicalSequence):
    """
    Class representing a protein (amino acid) sequence.

    Implements alphabet validation for 20 standard amino acids
    and provides methods for protein-specific analysis.
    """

    _valid_amino_acids: set[str] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*'
    }

    _molecular_weights: dict[str, float] = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
        'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.18,
        'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
        '*': 0.0
    }

    def __init__(self, sequence: str):
        """
        Initialize an amino acid sequence.

        Arguments:
            sequence: The raw amino acid sequence string.

        Raises:
            ValueError: If sequence contains invalid amino acid characters.
        """
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError(f"Invalid characters for AminoAcidSequence: {self._sequence}")

    def check_alphabet(self) -> bool:
        """Validate that sequence contains only valid amino acid characters."""
        return set(self._sequence).issubset(self._valid_amino_acids)

    def calculate_molecular_weight(self) -> float:
        """
        Calculate the molecular weight of the protein sequence.

        Returns:
            Molecular weight in Daltons (g/mol).
        """
        return sum(
            self._molecular_weights.get(aa, 0.0) for aa in self._sequence
        )


def normalize_bounds(bounds: Union[int, float, Tuple[float, float]]) -> Tuple[float, float]:
    """
    Normalize bounds to a (min, max) tuple.

    Arguments:
        bounds: Single number (interpreted as upper bound) or pair (min, max).

    Returns:
        Tuple of two floats (min, max).

    Raises:
        ValueError: If input is invalid.
    """
    if isinstance(bounds, (int, float)):
        return 0.0, float(bounds)
    elif isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        return float(bounds[0]), float(bounds[1])
    else:
        raise ValueError("Bounds must be a number or a pair of numbers")


def _prepare_output_path(output_fastq: str) -> str:
    """
    Prepare output path by creating directory and checking for existing files.

    Arguments:
    output_fastq: output filename

    Returns:
    Full output path

    Raises:
        FileExistsError: If output path already exists
    """
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_fastq)

    if os.path.exists(output_path):
        if os.path.isfile(output_path):
            raise FileExistsError(f"Output file already exists: {output_path}")
        else:
            raise FileExistsError(f"Path already exists and is not a file: {output_path}")

    return output_path


def filter_fastq(
        input_fastq: str,
        output_fastq: str = "filtered.fastq",
        gc_bounds: Union[int, float, Tuple[float, float]] = (0, 100),
        length_bounds: Union[int, float, Tuple[int, int]] = (0, 2 ** 32),
        quality_threshold: float = 0.0
) -> None:
    """
    Filter FASTQ reads using Biopython.

    Loads all records into memory, applies filtering, and writes results.

    Arguments:
        input_fastq: Path to input FASTQ file.
        output_fastq: Path to output filtered FASTQ file.
        gc_bounds: GC content bounds (min, max) as percentage.
        length_bounds: Sequence length bounds (min, max).
        quality_threshold: Minimum average Phred quality score.

    Raises:
        FileNotFoundError: If input file does not exist.
        FileExistsError: If output file already exists.
    """
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Input file not found: {input_fastq}")

    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    output_path = _prepare_output_path(output_fastq)

    filtered_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        seq_len = len(record.seq)

        if not (len_min <= seq_len <= len_max):
            continue

        gc_content = gc_fraction(record.seq) * 100
        if not (gc_min <= gc_content <= gc_max):
            continue

        mean_quality = sum(record.letter_annotations["phred_quality"]) / seq_len
        if mean_quality < quality_threshold:
            continue

        filtered_records.append(record)

    SeqIO.write(filtered_records, output_path, "fastq")


def filter_fastq_stream(
        input_fastq: str,
        output_fastq: str = "filtered.fastq",
        gc_bounds: Union[int, float, Tuple[float, float]] = (0, 100),
        length_bounds: Union[int, float, Tuple[int, int]] = (0, 2 ** 32),
        quality_threshold: float = 0.0
) -> None:
    """
    Filter FASTQ reads in streaming mode using Biopython.

    Processes one read at a time without loading all records into memory.
    Memory-efficient for large FASTQ files (>1 GB).

    Arguments:
        input_fastq: Path to input FASTQ file.
        output_fastq: Path to output filtered FASTQ file.
        gc_bounds: GC content bounds (min, max) as percentage.
        length_bounds: Sequence length bounds (min, max).
        quality_threshold: Minimum average Phred quality score.

    Raises:
        FileNotFoundError: If input file does not exist.
        FileExistsError: If output file already exists.
    """
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Input file not found: {input_fastq}")

    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    output_path = _prepare_output_path(output_fastq)

    def filtered_records():
        for record in SeqIO.parse(input_fastq, "fastq"):
            seq_len = len(record.seq)

            if not (len_min <= seq_len <= len_max):
                continue

            gc_content = gc_fraction(record.seq) * 100
            if not (gc_min <= gc_content <= gc_max):
                continue

            mean_quality = sum(record.letter_annotations["phred_quality"]) / seq_len
            if mean_quality < quality_threshold:
                continue

            yield record

    SeqIO.write(filtered_records(), output_path, "fastq")