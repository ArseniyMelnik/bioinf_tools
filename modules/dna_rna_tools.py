# Define nucleotide sets for DNA and RNA
DNA_BASES = {"A", "T", "G", "C"}
RNA_BASES = {"A", "U", "G", "C"}


def is_nucleic_acid(sequence):
    """
    Check if the given sequence is a valid DNA or RNA sequence.
    Returns True if all characters belong to DNA_BASES or RNA_BASES.
    """
    sequence_upper = sequence.upper()
    return (
        set(sequence_upper).issubset(DNA_BASES)
        or set(sequence_upper).issubset(RNA_BASES)
    )


def transcribe(sequence):
    """
    Transcribe a DNA sequence into RNA by replacing 'T' with 'U'.
    Case is preserved. Returns None if sequence is invalid.
    """
    if not is_nucleic_acid(sequence):
        print("Invalid sequence")
        return None

    # Replace 'T' with 'U' and 't' with 'u'; leave other characters unchanged
    return "".join(
        "U" if n == "T" else "u" if n == "t" else n for n in sequence
    )


def reverse(sequence):
    """
    Return the reversed sequence.
    Returns None if sequence is invalid.
    """
    if not is_nucleic_acid(sequence):
        print("Invalid sequence")
        return None
    return sequence[::-1]


def complement(sequence):
    """
    Return the complementary sequence.
    Supports both DNA and RNA, preserving case.
    Returns None if sequence is invalid.
    """
    if not is_nucleic_acid(sequence):
        print("Invalid sequence")
        return None

    # Dictionaries for DNA and RNA complements
    DNA_COMP = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }
    RNA_COMP = {
        "A": "U", "U": "A", "G": "C", "C": "G",
        "a": "u", "u": "a", "g": "c", "c": "g",
    }

    # Determine if sequence is RNA (contains 'U' or 'u') or DNA
    comp_dict = RNA_COMP if "U" in sequence or "u" in sequence else DNA_COMP

    # Return the complementary sequence
    return "".join(comp_dict[n] for n in sequence)


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA or RNA sequence.
    Combines the reverse and complement functions.
    Returns None if sequence is invalid.
    """
    if not is_nucleic_acid(sequence):
        print("Invalid sequence")
        return None
    return complement(reverse(sequence))


def run_dna_rna_tools(*args):
    """
    Universal wrapper function to apply a selected procedure
    on one or more sequences.

    Usage:
        run_dna_rna_tools(seq1, seq2, ..., procedure_name)
    Returns a single result if one sequence is provided,
    or a list of results for multiple sequences.
    """
    *sequences, procedure = args

    # Map procedure names to functions
    func_map = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    func = func_map[procedure]

    # Apply the selected function to all sequences
    results = [func(seq) for seq in sequences]

    # Return single result if only one sequence, else return list
    return results[0] if len(results) == 1 else results

