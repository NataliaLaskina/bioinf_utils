"""
Tests for bioinf_utils package.

Covers sequence classes, filtering functions, file operations,
and error handling.
"""

import pytest
import os
import tempfile
from bioinf_utils import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    filter_fastq,
    filter_fastq_stream,
    normalize_bounds
)
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    parse_blast_output
)


class TestDNASequence:
    """Tests for DNASequence class."""

    def test_dna_complement(self):
        """Test DNA complement calculation."""
        dna = DNASequence("ATGC")
        assert str(dna.complement()) == "TACG"

    def test_dna_reverse_complement(self):
        """Test DNA reverse complement calculation."""
        dna = DNASequence("ATGC")
        assert str(dna.reverse_complement()) == "GCAT"

    def test_dna_transcribe(self):
        """Test DNA to RNA transcription."""
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert str(rna) == "AUGC"
        assert isinstance(rna, RNASequence)

    def test_dna_slice_returns_object(self):
        """Test that slicing returns a DNASequence object."""
        dna = DNASequence("ATGCATGC")
        slice_result = dna[1:5]
        assert isinstance(slice_result, DNASequence)
        assert str(slice_result) == "TGCA"
        assert str(slice_result.complement()) == "ACGT"


class TestRNASequence:
    """Tests for RNASequence class."""

    def test_rna_complement(self):
        """Test RNA complement calculation."""
        rna = RNASequence("AUGC")
        assert str(rna.complement()) == "UACG"

    def test_rna_reverse_complement(self):
        """Test RNA reverse complement calculation."""
        rna = RNASequence("AUGC")
        assert str(rna.reverse_complement()) == "GCAU"


class TestAminoAcidSequence:
    """Tests for AminoAcidSequence class."""

    def test_molecular_weight(self):
        """Test protein molecular weight calculation."""
        protein = AminoAcidSequence("MKTAY")
        weight = protein.calculate_molecular_weight()
        assert isinstance(weight, float)
        assert weight > 0

    def test_amino_acid_valid_alphabet(self):
        """Test amino acid alphabet validation."""
        protein = AminoAcidSequence("MKTAY")
        assert protein.check_alphabet() is True


class TestErrorHandling:
    """Tests for error handling."""

    def test_invalid_dna_sequence_raises_value_error(self):
        """Test that invalid DNA sequence raises ValueError."""
        with pytest.raises(ValueError):
            DNASequence("ATGX")

    def test_file_not_found_raises_file_not_found_error(self):
        """Test that missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent.fastq", "output.fastq")


class TestFileOperations:
    """Tests for file operations."""

    def test_fasta_conversion_creates_oneline_file(self):
        """Test FASTA multiline to oneline conversion with file I/O."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGC\n")
            f.write("ATGC\n")
            input_file = f.name

        output_file = input_file.replace('.fasta', '_oneline.fasta')

        try:
            convert_multiline_fasta_to_oneline(input_file, output_file)

            assert os.path.exists(output_file)
            with open(output_file, 'r') as f:
                lines = f.readlines()
                assert len(lines) == 2
                assert lines[0].strip() == ">seq1"
                assert lines[1].strip() == "ATGCATGC"
        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_parse_blast_output_creates_sorted_file(self):
        """Test BLAST parsing creates sorted output file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("Query= test\n")
            f.write("Sequences producing significant alignments:\n")
            f.write(">MULTISPECIES: Protein A [Homo sapiens]\n")
            f.write(">Protein B [Mus musculus]\n")
            input_file = f.name

        output_file = input_file.replace('.txt', '_hits.txt')

        try:
            parse_blast_output(input_file, output_file)

            assert os.path.exists(output_file)
            with open(output_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
                assert lines == ["Protein A", "Protein B"]
        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_normalize_bounds_single_value(self):
        """Test bounds normalization with single value."""
        assert normalize_bounds(50) == (0.0, 50.0)

    def test_normalize_bounds_tuple(self):
        """Test bounds normalization with tuple."""
        assert normalize_bounds((20, 80)) == (20.0, 80.0)

    def test_normalize_bounds_list(self):
        """Test bounds normalization with list."""
        assert normalize_bounds([0, 100]) == (0.0, 100.0)