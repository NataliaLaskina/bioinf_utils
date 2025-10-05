"""
Main module for bioinf_utils package.
"""

from typing import Dict, Tuple, Union
from bioinf_utils.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)
from bioinf_utils.fastq_utils import (
    normalize_bounds,
    calculate_gc,
    calculate_mean_quality
)


def run_dna_rna_tools(*args: str) -> None:
    """
    Apply operation to one or more nucleotide sequences.

    Arguments:
    *args: sequences followed by operation name
           valid operations: 'is_nucleic_acid', 'transcribe',
           'reverse', 'complement', 'reverse_complement'

    Returns None. Prints result to stdout.
    Raises ValueError on invalid input.
    """
    if len(args) < 2:
        raise ValueError("Need at least one sequence and one operation")

    *sequences, operation = args

    ops = {
        'is_nucleic_acid': is_nucleic_acid,
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement
    }

    if operation not in ops:
        raise ValueError("Unknown operation")

    func = ops[operation]
    results = []

    for seq in sequences:
        if operation != 'is_nucleic_acid' and not is_nucleic_acid(seq):
            raise ValueError(f"Invalid nucleic acid sequence: {seq}")
        results.append(func(seq))

    if len(sequences) == 1:
        print(results[0])
    else:
        print(results)


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[int, float, Tuple[float, float]] = (0, 100),
    length_bounds: Union[int, float, Tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0.0
) -> Dict[str, Tuple[str, str]]:
    """
    Filter FASTQ reads by GC content, length, and quality.

    Arguments:
    seqs: dictionary {read_name: (sequence, quality_string)}
    gc_bounds: GC content bounds in percent (single value → upper bound)
    length_bounds: sequence length bounds (single value → upper bound)
    quality_threshold: minimum average Phred+33 quality

    Returns filtered dictionary of reads.
    """
    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    filtered: Dict[str, Tuple[str, str]] = {}

    for name, (seq, quality) in seqs.items():
        seq_len = len(seq)
        if not (len_min <= seq_len <= len_max):
            continue

        gc = calculate_gc(seq)
        if not (gc_min <= gc <= gc_max):
            continue

        mean_qual = calculate_mean_quality(quality)
        if mean_qual < quality_threshold:
            continue

        filtered[name] = (seq, quality)

    print(filtered)
