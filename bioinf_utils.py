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
