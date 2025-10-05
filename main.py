# main.py

from modules.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)

from modules.filter_fastq import filter_fastq

def run_fastq_tools(seqs, gc_bounds=(0,100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    Universal wrapper to call filter_fastq with optional filtering parameters.
    """
    return filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)


# run_dna_rna_tools
def run_dna_rna_tools(*args):
    """
    Universal function for working with sequences.
    Accepts any number of sequences followed by procedure name.
    """
    *sequences, procedure = args
    func = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }[procedure]

    results = [func(seq) for seq in sequences]
    return results[0] if len(results) == 1 else results


