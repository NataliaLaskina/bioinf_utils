def is_nucleic_acid(seq: str) -> bool:
    """
    Check if sequence contains only valid DNA or RNA nucleotides.

    Arguments:
    seq: input nucleotide sequence (case-insensitive)

    Returns bool.
    """
    if not seq:
        return False
    bases = set(seq.upper())
    valid_dna = {'A', 'T', 'G', 'C'}
    valid_rna = {'A', 'U', 'G', 'C'}
    return bases <= valid_dna or bases <= valid_rna


def transcribe(seq: str) -> str:
    """
    Transcribe DNA to RNA by replacing T with U.

    Arguments:
    seq: DNA sequence (may contain lowercase letters)

    Returns RNA sequence.
    """
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq: str) -> str:
    """
    Return the reverse of a sequence.

    Arguments:
    seq: input sequence

    Returns reversed sequence.
    """
    return seq[::-1]