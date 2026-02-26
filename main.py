from modules.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)
from modules.filter_fastq import filter_fastq

# Call filter_fastq directly
# DNA/RNA tools wrapper
def run_dna_rna_tools(*args):
    """
    Universal wrapper for DNA/RNA procedures.
    Accepts any number of sequences followed by procedure name.
    """
    *sequences, procedure = args

    func_map = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement
    }

    func = func_map[procedure]
    results = [func(seq) for seq in sequences]
    return results[0] if len(results) == 1 else results
