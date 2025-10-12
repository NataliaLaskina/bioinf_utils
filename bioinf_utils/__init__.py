"""
Bioinformatics utilities package.
Provides tools for DNA/RNA sequence manipulation and FASTQ file processing,
including reading, writing, and quality/GC/length analysis.
"""

from .dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)

from .fastq_utils import (
    normalize_bounds,
    calculate_gc,
    calculate_mean_quality,
    is_gc_within_bounds,
    is_length_within_bounds,
    is_quality_above_threshold,
    read_fastq,
    write_fastq
)

__all__ = [
    "is_nucleic_acid",
    "transcribe",
    "reverse",
    "complement",
    "reverse_complement",
    "normalize_bounds",
    "calculate_gc",
    "calculate_mean_quality",
    "is_gc_within_bounds",
    "is_length_within_bounds",
    "is_quality_above_threshold",
    "read_fastq",
    "write_fastq"
]