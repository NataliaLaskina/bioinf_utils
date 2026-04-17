#!/usr/bin/env python3
"""
Command-line interface for bioinf_utils package.

Provides CLI access to sequence manipulation, FASTQ filtering,
FASTA conversion, and BLAST output parsing.
"""

import argparse
import sys
import logging

from bioinf_utils import (
    filter_fastq,
    filter_fastq_stream,
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    normalize_bounds
)
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    parse_blast_output
)

# Configure logging for CLI usage
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('bioinf_utils.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments using argparse."""
    parser = argparse.ArgumentParser(
        description="Bioinformatics utilities for sequence manipulation and file processing",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # === filter-fastq command ===
    filter_parser = subparsers.add_parser(
        "filter-fastq",
        help="Filter FASTQ file by GC content, length, and quality"
    )
    filter_parser.add_argument("input", help="Path to input FASTQ file")
    filter_parser.add_argument("output", help="Path to output filtered FASTQ file")
    filter_parser.add_argument("--gc-min", type=float, default=0, help="Minimum GC content (%%)")
    filter_parser.add_argument("--gc-max", type=float, default=100, help="Maximum GC content (%%)")
    filter_parser.add_argument("--len-min", type=int, default=0, help="Minimum sequence length")
    filter_parser.add_argument("--len-max", type=int, default=2 ** 32, help="Maximum sequence length")
    filter_parser.add_argument("--quality", type=float, default=0.0, help="Minimum average Phred quality")
    filter_parser.add_argument("--stream", action="store_true", help="Use memory-efficient streaming mode")

    # === convert-fasta command ===
    fasta_parser = subparsers.add_parser(
        "convert-fasta",
        help="Convert multiline FASTA to oneline format"
    )
    fasta_parser.add_argument("input", help="Path to input FASTA file")
    fasta_parser.add_argument("--output", help="Path to output FASTA file (default: <input>_oneline.fasta)")

    # === parse-blast command ===
    blast_parser = subparsers.add_parser(
        "parse-blast",
        help="Parse BLAST output and extract protein names"
    )
    blast_parser.add_argument("input", help="Path to BLAST results file")
    blast_parser.add_argument("output", help="Path to output file with protein names")

    # === seq-info command ===
    seq_parser = subparsers.add_parser(
        "seq-info",
        help="Analyze biological sequence"
    )
    seq_parser.add_argument("sequence", help="Sequence string")
    seq_parser.add_argument(
        "--type",
        choices=["dna", "rna", "protein"],
        required=True,
        help="Type of biological sequence"
    )

    return parser.parse_args()


def handle_filter_fastq(args: argparse.Namespace) -> None:
    """Execute FASTQ filtering based on CLI arguments."""
    gc_bounds = (args.gc_min, args.gc_max)
    length_bounds = (args.len_min, args.len_max)

    logger.info(f"Filtering FASTQ: {args.input} -> {args.output}")
    logger.info(f"Parameters: gc={gc_bounds}, len={length_bounds}, quality>={args.quality}, stream={args.stream}")

    if args.stream:
        filter_fastq_stream(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=args.quality
        )
    else:
        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=args.quality
        )
    logger.info("FASTQ filtering completed successfully")


def handle_convert_fasta(args: argparse.Namespace) -> None:
    """Execute FASTA conversion based on CLI arguments."""
    logger.info(f"Converting FASTA: {args.input} -> {args.output}")
    convert_multiline_fasta_to_oneline(
        input_fasta=args.input,
        output_fasta=args.output
    )
    logger.info("FASTA conversion completed successfully")


def handle_parse_blast(args: argparse.Namespace) -> None:
    """Execute BLAST parsing based on CLI arguments."""
    logger.info(f"Parsing BLAST: {args.input} -> {args.output}")
    parse_blast_output(
        input_file=args.input,
        output_file=args.output
    )
    logger.info("BLAST parsing completed successfully")


def handle_seq_info(args: argparse.Namespace) -> None:
    """Execute sequence analysis based on CLI arguments."""
    seq_str = args.sequence.upper()
    logger.info(f"Analyzing {args.type} sequence: {seq_str}")

    if args.type == "dna":
        seq = DNASequence(seq_str)
        print(f"Type: DNA")
        print(f"Length: {len(seq)}")
        print(f"Complement: {seq.complement()}")
        print(f"Reverse complement: {seq.reverse_complement()}")
        print(f"Transcribed RNA: {seq.transcribe()}")
        print(f"Valid alphabet: {seq.check_alphabet()}")

    elif args.type == "rna":
        seq = RNASequence(seq_str)
        print(f"Type: RNA")
        print(f"Length: {len(seq)}")
        print(f"Complement: {seq.complement()}")
        print(f"Reverse complement: {seq.reverse_complement()}")
        print(f"Valid alphabet: {seq.check_alphabet()}")

    elif args.type == "protein":
        seq = AminoAcidSequence(seq_str)
        print(f"Type: Protein")
        print(f"Length: {len(seq)}")
        print(f"Molecular weight: {seq.calculate_molecular_weight():.2f} Da")
        print(f"Valid alphabet: {seq.check_alphabet()}")

    logger.info("Sequence analysis completed successfully")


def main() -> None:
    """Main entry point for CLI."""
    args = parse_args()

    if args.command is None:
        print("Error: No command specified. Use -h for help.")
        sys.exit(1)

    try:
        if args.command == "filter-fastq":
            handle_filter_fastq(args)
        elif args.command == "convert-fasta":
            handle_convert_fasta(args)
        elif args.command == "parse-blast":
            handle_parse_blast(args)
        elif args.command == "seq-info":
            handle_seq_info(args)

        print("\n✓ Operation completed successfully!")

    except FileNotFoundError as e:
        logger.error(str(e))
        print(f"Error: {e}")
        sys.exit(1)
    except FileExistsError as e:
        logger.error(str(e))
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(str(e))
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()