# Bioinformatics Sequence Tools

## Classes

### `BiologicalSequence`

* Abstract base class for working with biological sequences.
* Supports:

  * `len(seq)` — sequence length
  * `seq[index]` and `seq[start:stop]` — indexing and slicing
  * `print(seq)` — nicely formatted output
  * `check_alphabet()` — validates sequence characters

### `NucleicAcidSequence`

* Parent class for DNA and RNA sequences.
* Methods:

  * `complement()` — returns complementary sequence
  * `reverse()` — returns reversed sequence
  * `reverse_complement()` — returns reverse-complemented sequence
* **Direct instantiation is not recommended.** Use subclasses `DNASequence` or `RNASequence`.

### `DNASequence` (inherits from NucleicAcidSequence)

* Additional method:

  * `transcribe()` — returns an `RNASequence` object

### `RNASequence` (inherits from NucleicAcidSequence)

* Inherits complement and reverse methods from the parent class

### `AminoAcidSequence`

* Methods:

  * `aa_composition()` — returns the percentage composition of each amino acid
  * `len(seq)` — sequence length
  * `seq[start:stop]` — slicing
  * Additional methods can include hydrophobicity, molecular weight, etc.

---

## Usage Examples

```python
from main import DNASequence, RNASequence, AminoAcidSequence

# DNA example
dna = DNASequence("ATGCGT")
print("Complement:", dna.complement())
rna = dna.transcribe()
print("Reverse complement of RNA:", rna.reverse_complement())

# Amino acid example
protein = AminoAcidSequence("ACDEFGHIK")
print("Amino acid composition (%):", protein.aa_composition())
print("Length:", len(protein))
print("Slice (0:5):", protein[:5])
```

---

## FASTQ Filtering Function

### `filter_fastq(input_fastq, output_fastq, gc_bounds=(0,100), length_bounds=(0,inf), quality_threshold=0)`

Filters FASTQ files using **Biopython** according to:

* **GC content**: Only sequences with GC percentage within `gc_bounds` are kept.
* **Sequence length**: Only sequences with length within `length_bounds` are kept.
* **Average Phred quality**: Only sequences with mean quality ≥ `quality_threshold` are kept.

**Parameters:**

* `input_fastq` — path to the input FASTQ file (e.g., `data/example.fastq`)
* `output_fastq` — path for the filtered FASTQ file (e.g., `filtered/example_filtered.fastq`)
* `gc_bounds` — tuple `(min, max)` or single value for maximum GC %
* `length_bounds` — tuple `(min, max)` or single value for maximum sequence length
* `quality_threshold` — minimum average Phred quality score

**Returns:**

A dictionary with statistics and filtered records:

```python
{
    "total": total_reads,           # total reads processed
    "passed": passed_reads,         # reads that passed all filters
    "filtered_records": {           # dictionary of SeqRecord objects that passed
        sequence_id: SeqRecord,
        ...
    }
}
```

**Example:**

```python
from main import filter_fastq

result = filter_fastq(
    input_fastq="data/example.fastq",
    output_fastq="filtered/example_filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 150),
    quality_threshold=20
)

print("Total reads:", result["total"])
print("Reads passed filter:", result["passed"])
```

**Notes:**

* Input FASTQ files should have matching sequence and quality string lengths.
* Output FASTQ files are automatically created in the specified folder (`filtered/` in this repository).
* Requires **Biopython**: `pip install biopython` or 'pip install -r requirements.txt'

---

