# bioinf_tools

**bioinf_tools** is a set of tools for working with nucleic acids, FASTQ files, and bioinformatics data. The package provides functionality for sequence analysis, filtering, and file processing.
---

### 1. DNA/RNA Sequence Operations (`dna_rna_tools.py`)
- Check if a sequence is valid DNA or RNA
- Transcribe DNA to RNA (preserving case)
- Reverse a sequence
- Generate complementary sequences
- Generate reverse-complement sequences

### 2. FASTQ Filtering (`filter_fastq.py`)
- Filter sequences by GC content
- Filter sequences by length
- Filter sequences by mean Phred33 quality
- Read FASTQ files into a dictionary `{name: (sequence, quality)}`
- Write filtered sequences back to FASTQ files
- Safe file handling (prevents overwriting existing files)

### 3. Bioinformatics File Processing (`bio_files_processor.py`)
- Convert multi-line FASTA files into single-line FASTA
- Parse BLAST txt output to extract top hit descriptions
---

## Installation

Clone the repository:

```bash
git clone https://github.com/ArseniyMelnik/bioinf_tools.git
cd bioinf_tools
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
from main import filter_fastq, read_fastq, write_fastq

# Read sequences from a FASTQ file
seqs = read_fastq("data/example.fastq")

# Filter sequences by GC content, length, and quality
filtered = filter_fastq(seqs, gc_bounds=(40,60), length_bounds=(5,1000), quality_threshold=30)

# Save filtered sequences to a new FASTQ file
write_fastq(filtered, "filtered/filtered_output.fastq")

# Parameters for filtering:
# - gc_bounds — tuple (min, max) or single number (upper bound)
# - length_bounds — tuple (min, max) or single number (upper bound)
# - quality_threshold — minimum mean Phred33 quality
```

### Bioinformatics File Processing

```python
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output

# Convert multiline FASTA to single-line FASTA
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")

# Parse BLAST output to extract top hits
parse_blast_output("blast_results.txt", "top_hits.txt")
```
### Notes on Data Folders

data/ — this directory is used to store input biological data files, such as FASTQ or FASTA
You can place your test or working datasets here.

filtered/ — this folder is automatically created when running filtering functions (e.g., filter_fastq).
It contains the output files produced after sequence filtering.
The script ensures that:

existing files will not be overwritten accidentally;

the folder is created if it doesn’t exist.

## Project Structure
```
bioinf_tools/
│
├── main.py                # Main interface script
├── bio_files_processor.py # Additional bioinformatics file processing utilities
├── README.md              # Project documentation
│
├── modules/               # Internal modules
│   ├── dna_rna_tools.py   # DNA/RNA sequence operations
│   └── filter_fastq.py    # FASTQ filtering utilities
│
├── data/                  # Input files (e.g., .fastq, .fasta, .gbk)
│   └── example.fastq
│
└── filtered/              # Automatically created folder for filtered outputs
    └── example_filtered.fastq

```

## Contact

Author: Arseniy Melnik
Email: 14amelnikd@gmail.com

GitHub: ArseniyMelnik
