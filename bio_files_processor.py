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