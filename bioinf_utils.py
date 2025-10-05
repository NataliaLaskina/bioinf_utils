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
