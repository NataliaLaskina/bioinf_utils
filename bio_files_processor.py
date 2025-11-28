import os
from typing import Optional


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: Optional[str] = None
) -> None:
    """
    Convert a multiline FASTA file to oneline format.

    Arguments:
    input_fasta: path to input FASTA file (multiline allowed)
    output_fasta: optional path to output file.
                  If not provided, uses '<input>_oneline.fasta'

    Reads sequences that may span multiple lines and writes each as a single line.
    """
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

    if output_fasta is None:
        base_name = os.path.splitext(os.path.basename(input_fasta))[0]
        output_fasta = f"{base_name}_oneline.fasta"

    if os.path.exists(output_fasta):
        if os.path.isfile(output_fasta):
            raise FileExistsError(f"Output file already exists: {output_fasta}")
        else:
            raise FileExistsError(f"Path exists and is not a file: {output_fasta}")

    current_header = None
    current_seq_lines = []

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_header is not None:
                    full_seq = ''.join(current_seq_lines)
                    outfile.write(f"{current_header}\n{full_seq}\n")
                current_header = line
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        if current_header is not None:
            full_seq = ''.join(current_seq_lines)
            outfile.write(f"{current_header}\n{full_seq}\n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parse BLAST output file and extract top hit descriptions.

    Arguments:
    input_file: path to BLAST results in text format
    output_file: path to output file with sorted protein names (one per line)

    Extracts the first hit description for each query,
    removes species in brackets and common prefixes (e.g., 'MULTISPECIES: '),
    and saves all results (including duplicates) sorted alphabetically.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input BLAST file not found: {input_file}")

    protein_names = []  # Preserve duplicates
    prefixes_to_remove = ["MULTISPECIES: ", "PREDICTED: "]

    in_significant_alignments = False

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if "Sequences producing significant alignments:" in line:
                in_significant_alignments = True
                continue

            if in_significant_alignments and line.startswith('>'):
                desc = line[1:].strip()

                # Remove known non-protein-name prefixes
                for prefix in prefixes_to_remove:
                    if desc.startswith(prefix):
                        desc = desc[len(prefix):]
                        break  # only one prefix expected

                # Remove species in brackets, e.g. [Escherichia coli]
                if '[' in desc:
                    protein_name = desc.split('[', 1)[0].strip()
                else:
                    protein_name = desc.strip()

                protein_name = protein_name.rstrip('.').strip()
                if protein_name:
                    protein_names.append(protein_name)

                in_significant_alignments = False

    # Sort all names (duplicates preserved)
    sorted_proteins = sorted(protein_names)

    with open(output_file, 'w') as out_f:
        for name in sorted_proteins:
            out_f.write(f"{name}\n")