from typing import Union, Tuple


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
