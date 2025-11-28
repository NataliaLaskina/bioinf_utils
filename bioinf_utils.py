import os
from typing import Union, Tuple
from bioinf_utils.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)
from bioinf_utils.fastq_utils import (
    read_fastq,
    write_fastq,
    normalize_bounds,
    is_gc_within_bounds,
    is_length_within_bounds,
    is_quality_above_threshold
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


def _prepare_output_path(output_fastq: str) -> str:
    """
    Prepare output path by creating directory and checking for existing files.

    Arguments:
    output_fastq: output filename

    Returns:
    Full output path
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
    Filter FASTQ reads from a file and save results to another file.
    Loads entire file into memory.
    """
    seqs = read_fastq(input_fastq)

    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    filtered = {}
    for name, (seq, plus_line, quality) in seqs.items():
        if not is_length_within_bounds(seq, (len_min, len_max)):
            continue
        if not is_gc_within_bounds(seq, (gc_min, gc_max)):
            continue
        if not is_quality_above_threshold(quality, quality_threshold):
            continue
        filtered[name] = (seq, plus_line, quality)

    output_path = _prepare_output_path(output_fastq)
    write_fastq(filtered, output_path)


def filter_fastq_stream(
        input_fastq: str,
        output_fastq: str = "filtered.fastq",
        gc_bounds: Union[int, float, Tuple[float, float]] = (0, 100),
        length_bounds: Union[int, float, Tuple[int, int]] = (0, 2 ** 32),
        quality_threshold: float = 0.0
) -> None:
    """
    Filter FASTQ reads in streaming mode (memory-efficient).
    Processes one read at a time without loading all into memory.
    """
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Input file not found: {input_fastq}")

    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    output_path = _prepare_output_path(output_fastq)

    with open(input_fastq, 'r') as infile, open(output_path, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: expected '@', got {header}")
            seq = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            if not (seq and plus.startswith('+') and quality):
                raise ValueError("Invalid FASTQ format: incomplete record")

            if not is_length_within_bounds(seq, (len_min, len_max)):
                continue
            if not is_gc_within_bounds(seq, (gc_min, gc_max)):
                continue
            if not is_quality_above_threshold(quality, quality_threshold):
                continue

            outfile.write(f"{header}\n{seq}\n{plus}\n{quality}\n")