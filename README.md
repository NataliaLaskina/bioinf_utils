# bioinf_utils

A lightweight Python package for basic bioinformatics operations:  
DNA/RNA sequence manipulation and FASTQ read filtering.

## Installation

Clone the repository and use the modules directly (no external dependencies):

```bash
git clone https://github.com/NataliaLaskina/bioinf_utils.git
cd bioinf_utils
```
> ✅ Works with Python 3.9+

## Usage

### DNA/RNA Tools

Use run_dna_rna_tools() to apply operations to sequences:

```Python
from bioinf_utils import run_dna_rna_tools

# Transcribe DNA to RNA
run_dna_rna_tools("ATGC", "transcribe")  # Output: AUG C

# Get reverse complement
run_dna_rna_tools("ATGC", "reverse_complement")  # Output: GCAT
```

**Supported operations**:  
`transcribe`, `reverse`, `complement`, `reverse_complement`, `is_nucleic_acid`

### FASTQ Filtering

Filter reads by GC content, length, and quality:

```Python
from bioinf_utils import filter_fastq

seqs = {
    "read1": ("ATGC", "IIII"),
    "read2": ("AAAA", "!!!!")
}

filtered = filter_fastq(
    seqs,
    gc_bounds=(20, 80),      # GC% between 20 and 80
    length_bounds=(3, 10),   # Length between 3 and 10
    quality_threshold=30     # Mean Phred+33 quality ≥ 30
)
```

## Functions

- `run_dna_rna_tools(*sequences, operation)`  
  Applies a transformation to one or more nucleotide sequences.

- `filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)`  
  Filters FASTQ-like data based on biological and quality criteria.


## Author

Natalia Laskina 
lask.natalia@gmail.com  
