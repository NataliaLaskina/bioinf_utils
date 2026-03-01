# bioinf_utils

A lightweight Python package for common bioinformatics tasks: DNA/RNA/protein sequence manipulation (OOP), FASTQ filtering (Biopython), FASTA/BLAST parsing, and more.  
Designed for educational purposes and small-scale analysis.

> ✅ Works with Python 3.9+  
> ✅ Requires `biopython>=1.81`

## Installation

Clone the repository and install dependencies:

```
bash
git clone https://github.com/NataliaLaskina/bioinf_utils.git  
cd bioinf_utils
pip install -r requirements.txt
```

All functionality is available by importing classes and functions into your own scripts.

## Usage

### Object-Oriented Sequence Manipulation

#### Classes Overview

| Class | Description | Key Methods |
|-------|-------------|-------------|
| `DNASequence` | DNA nucleotide sequence | `complement()`, `reverse()`, `reverse_complement()`, `transcribe()` |
| `RNASequence` | RNA nucleotide sequence | `complement()`, `reverse()`, `reverse_complement()` |
| `AminoAcidSequence` | Protein (amino acid) sequence | `calculate_molecular_weight()`, `check_alphabet()` |

All classes inherit from `BiologicalSequence` and support:
- `len(seq)` — sequence length
- `seq[i]`, `seq[i:j]` — indexing and slicing
- `str(seq)` — human-readable string representation
- `check_alphabet()` — validate sequence characters

#### DNASequence Example

```
python
from bioinf_utils import DNASequence

dna = DNASequence("ATGC")

print(dna.complement())           # TACG (DNASequence)
print(dna.reverse())              # CGTA (DNASequence)
print(dna.reverse_complement())   # GCAT (DNASequence)
print(dna.transcribe())           # AUGC (RNASequence)
print(dna[1:3])                   # TG
print(dna.check_alphabet())       # True
```

#### RNASequence Example

```
python
from bioinf_utils import RNASequence

rna = RNASequence("AUGC")

print(rna.complement())           # UACG (RNASequence)
print(rna.reverse_complement())   # GCAU (RNASequence)
```

#### AminoAcidSequence Example

```
python
from bioinf_utils import AminoAcidSequence

protein = AminoAcidSequence("MKTAY")

print(len(protein))                           # 5
print(protein.calculate_molecular_weight())   # 684.8 Da
print(protein.check_alphabet())               # True
```

### FASTQ Filtering (Biopython)

#### filter_fastq(input_fastq, output_fastq, gc_bounds=(0,100), length_bounds=(0,2**32), quality_threshold=0.0)

Filters reads from a FASTQ file based on:

- **GC content** (`gc_bounds`): e.g., `(20, 80)` or `44.4` (upper bound only)
- **Read length** (`length_bounds`): e.g., `(50, 300)`
- **Average Phred+33 quality** (`quality_threshold`)

✅ Loads all records into memory, applies filtering, and writes results to `filtered/{output_fastq}`.

**Example**:

```
python
from bioinf_utils import filter_fastq

filter_fastq(
    input_fastq="reads.fastq",
    output_fastq="filtered_reads.fastq",
    gc_bounds=(30, 70),
    length_bounds=(100, 500),
    quality_threshold=20.0
)
```

#### filter_fastq_stream(input_fastq, output_fastq, ...)

Memory-efficient streaming version: processes one read at a time without loading the whole file into RAM.

✅ Ideal for large FASTQ files (>1 GB).

⚠️ **Important**: Requires standard FASTQ formatting (4 lines per record, newline at end of each line).

**Example**:

```
python
from bioinf_utils import filter_fastq_stream

filter_fastq_stream(
    input_fastq="large_reads.fastq",
    output_fastq="filtered_reads.fastq",
    gc_bounds=(30, 70),
    length_bounds=(100, 500),
    quality_threshold=20.0
)
```

### FASTA Processing

#### convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None)

Converts a multiline FASTA file (where sequences span multiple lines) into oneline format (one sequence per line).

- If `output_fasta` is not provided, output is saved as `{input}_oneline.fasta`.
- Ensures no accidental overwriting of existing files.

**Example**:

```
python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("genome.fasta", "genome_oneline.fasta")
``

### BLAST Output Parsing

#### parse_blast_output(input_file, output_file)

Parses a BLAST result file (text format) and extracts the top hit description for each query.

- Extracts protein names from lines starting with `>`
- Removes species annotations (text in `[...]`)
- Removes common prefixes (`MULTISPECIES:`, `PREDICTED:`)
- Preserves all hits (one per query), including duplicates
- Outputs a sorted list of protein names (one per line)

**Ideal for**: Preparing a clean list of candidate proteins for downstream analysis.

**Example**:

```
python
from bio_files_processor import parse_blast_output

parse_blast_output("blast_results.txt", "top_hits.txt")
```

## Project Structure

```
bioinf_utils/
├── bioinf_utils.py              # Main module: OOP classes + FASTQ filtering
├── bio_files_processor.py       # FASTA conversion + BLAST parsing
├── requirements.txt             # Dependencies (biopython)
└── README.md                    # This file
```

## Dependencies

- Python 3.9+
- Biopython >= 1.81

Install with:
```
bash
pip install -r requirements.txt
```

## Author

Natalia Laskina  
lask.natalia@gmail.com