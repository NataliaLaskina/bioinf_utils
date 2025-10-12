import os
from typing import Dict, Tuple, Union, List, Optional


def normalize_bounds(bounds: Union[int, float, Tuple[float, float]]) -> Tuple[float, float]:
    """
    Normalize bounds to a (min, max) tuple.

    Arguments:
    bounds: single number (interpreted as upper bound) or pair (min, max)

    Returns tuple of two floats (min, max).
    Raises ValueError if input is invalid.
    """
    if isinstance(bounds, (int, float)):
        return 0.0, float(bounds)
    elif isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        return float(bounds[0]), float(bounds[1])
    else:
        raise ValueError("Bounds must be a number or a pair of numbers")


def calculate_gc(seq: str) -> float:
    """
    Calculate GC content percentage.

    Arguments:
    seq: nucleotide sequence (case-insensitive)

    Returns GC content as percentage (0.0 to 100.0).
    """
    if not seq:
        return 0.0
    gc_count = sum(1 for base in seq.upper() if base in ('G', 'C'))
    return (gc_count / len(seq)) * 100.0


def calculate_mean_quality(quality_str: str) -> float:
    """
    Calculate average Phred+33 quality score.

    Arguments:
    quality_str: FASTQ quality string (ASCII-encoded)

    Returns mean quality score.
    """
    if not quality_str:
        return 0.0
    total = sum(ord(char) - 33 for char in quality_str)
    return total / len(quality_str)


def is_gc_within_bounds(seq: str, gc_bounds: Tuple[float, float]) -> bool:
    """
    Check if GC content of sequence is within given bounds.

    Arguments:
    seq: nucleotide sequence
    gc_bounds: (min, max) GC percentage

    Returns bool.
    """
    gc = calculate_gc(seq)
    gc_min, gc_max = gc_bounds
    return gc_min <= gc <= gc_max


def is_length_within_bounds(seq: str, length_bounds: Tuple[float, float]) -> bool:
    """
    Check if sequence length is within given bounds.

    Arguments:
    seq: nucleotide sequence
    length_bounds: (min, max) length

    Returns bool.
    """
    l = len(seq)
    len_min, len_max = length_bounds
    return len_min <= l <= len_max


def is_quality_above_threshold(quality_str: str, threshold: float) -> bool:
    """
    Check if average quality is at or above threshold.

    Arguments:
    quality_str: FASTQ quality string
    threshold: minimum acceptable average Phred+33 score

    Returns bool.
    """
    mean_qual = calculate_mean_quality(quality_str)
    return mean_qual >= threshold


def read_fastq(file_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Read a FASTQ file and return its contents as a dictionary.

    Arguments:
    file_path: path to the input FASTQ file

    Returns:
    Dictionary {read_name: (sequence, quality_string)}
    Raises FileNotFoundError if file does not exist.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    seqs = {}
    with open(file_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: expected '@' at start of header, got {header}")
            seq = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            if not (seq and plus == '+' and quality):
                raise ValueError("Invalid FASTQ format: missing sequence or quality lines")

            read_name = header[1:]  # remove '@'
            seqs[read_name] = (seq, quality)

    return seqs


def write_fastq(seqs: Dict[str, Tuple[str, str]], output_filename: str) -> None:
    """
    Write filtered FASTQ sequences to a file in the 'filtered/' directory.

    Arguments:
    seqs: dictionary of filtered reads
    output_filename: name of the output file (without path)

    Creates 'filtered/' directory if it doesn't exist.
    Raises FileExistsError if output path already exists (as file or directory).
    """
    output_dir = "filtered"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)

    if os.path.exists(output_path):
        if os.path.isfile(output_path):
            raise FileExistsError(f"Output file already exists: {output_path}")
        else:
            raise FileExistsError(f"Path already exists and is not a file: {output_path}")

    with open(output_path, 'w') as f:
        for name, (seq, quality) in seqs.items():
            f.write(f"@{name}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{quality}\n")
            