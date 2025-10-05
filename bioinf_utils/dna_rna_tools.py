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
