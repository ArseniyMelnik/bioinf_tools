# bioinf_tools

This package provides tools for working with DNA/RNA sequences and FASTQ files. It allows filtering sequences by GC content, length, and quality, and performing basic operations on nucleic acids such as transcription, reverse, complement, and reverse complement.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/your-repo.git
cd your-repo
```
## Usage examples
### DNA/RNA tools
```python
from main import run_dna_rna_tools

seq = "ATGCGT"
result = run_dna_rna_tools(seq, "reverse_complement")
print(result)

# Available procedures:
# - is_nucleic_acid
# - transcribe
# - reverse
# - complement
# - reverse_complement
```
### FASTQ filtering

```python
from main import run_fastq_tools

seqs = {
    "read1": ("ATGCGT", "IIIIII"),
    "read2": ("GGCATG", "HHHHHH")
}

filtered = run_fastq_tools(seqs, gc_bounds=(40,60), quality_threshold=30)
print(filtered)

# Parameters for filtering:
# - gc_bounds — tuple (min, max) or single number (upper bound)
# - length_bounds — tuple (min, max) or single number (upper bound)
# - quality_threshold — minimum mean Phred33 quality
```

## Project Structure
your-repo/
│
├── main.py             # Main interface
├── modules/
│   ├── dna_rna_tools.py  # DNA/RNA sequence operations
│   └── filter_fastq.py   # FASTQ filtering utilities
└── README.md


## Contact

Author: Arseniy Melnik
Email: 14amelnikd@gmail.com

GitHub: ArseniyMelnik
