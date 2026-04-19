import abc
import os
import argparse
import logging
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# ======================================
# Logging setup
# ======================================
logging.basicConfig(
    filename="fastq_filter.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# ======================================
# BiologicalSequence
# ======================================
class BiologicalSequence(abc.ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence
        if not self.check_alphabet():
            raise ValueError(f"Sequence contains invalid characters: {sequence}")

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        # Index and slice support
        return self.__class__(self.sequence[key]) if isinstance(key, slice) else self.sequence[key]

    def __str__(self):
        return f"{self.__class__.__name__}: {self.sequence}"

    @abc.abstractmethod
    def check_alphabet(self) -> bool:
        """
        Checks whether the sequence contains only valid characters from the defined alphabet
        """
        pass

# ======================================
# NucleicAcidSequence
# ======================================
class NucleicAcidSequence(BiologicalSequence):
    # Complementary dictionaries for DNA and RNA
    _DNA_COMP = {"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g"}
    _RNA_COMP = {"A":"U","U":"A","G":"C","C":"G","a":"u","u":"a","g":"c","c":"g"}

    @abc.abstractmethod
    def _complement_table(self) -> dict:
        """
        Returns the complement mapping table for the specific type of nucleic acid sequence
        """
        raise NotImplementedError

    def complement(self):
        table = self._complement_table()
        return self.__class__("".join(table[n] for n in self.sequence))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.reverse().complement()

# ======================================
# DNASequence & RNASequence
# ======================================
class DNASequence(NucleicAcidSequence):
    _alphabet = {"A","T","G","C","a","t","g","c"}

    def check_alphabet(self) -> bool:
        return set(self.sequence).issubset(self._alphabet)

    def _complement_table(self):
        return self._DNA_COMP

    def transcribe(self):
        # Return the RNASequence object
        rna_seq = "".join("U" if n=="T" else "u" if n=="t" else n for n in self.sequence)
        return RNASequence(rna_seq)

class RNASequence(NucleicAcidSequence):
    _alphabet = {"A","U","G","C","a","u","g","c"}

    def check_alphabet(self) -> bool:
        return set(self.sequence).issubset(self._alphabet)

    def _complement_table(self):
        return self._RNA_COMP

# ======================================
# AminoAcidSequence
# ======================================
class AminoAcidSequence(BiologicalSequence):
    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")  # 20 стандартных аминокислот

    def check_alphabet(self) -> bool:
        return set(self.sequence.upper()).issubset(self._alphabet)

    def aa_composition(self):
        composition = {}
        for aa in self.sequence.upper():
            composition[aa] = composition.get(aa, 0) + 1
        total = len(self.sequence)
        return {aa: count/total*100 for aa, count in composition.items()}


# ======================================
# filter_fastq
# ======================================

def filter_fastq(input_fastq: str,
                 output_fastq: str,
                 gc_bounds=(0, 100),
                 length_bounds=(0, float("inf")),
                 quality_threshold=0) -> dict:
    """
    Filter a FastQ file using Biopython.
    
    Parameters:
        input_fastq: path to input FastQ file
        output_fastq: path to output FastQ file
        gc_bounds: tuple(min, max) or single int (max)
        length_bounds: tuple(min, max) or single int (max)
        quality_threshold: minimum average Phred quality score
    
    Returns:
        dict: {
            "total": total reads,
            "passed": reads passed filters,
            "filtered_records": {sequence_id: SeqRecord}
        }
    """

    logging.info(f"Started filtering: {input_fastq}")

    # input check
    if not os.path.exists(input_fastq):
        logging.error(f"File not found: {input_fastq}")
        raise FileNotFoundError(f"Input FASTQ file not found: {input_fastq}")

    # Normalize bounds
    gc_min, gc_max = (0, gc_bounds) if isinstance(gc_bounds, (int, float)) else gc_bounds
    len_min, len_max = (0, length_bounds) if isinstance(length_bounds, (int, float)) else length_bounds

    out_dir = os.path.dirname(output_fastq)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    filtered = {}

    # Counts
    total_count = 0
    passed_count = 0

    with open(output_fastq, "w") as out_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            total_count += 1

            seq_len = len(record.seq)
            if seq_len == 0:
                continue

            seq_gc = gc_fraction(record.seq) * 100
            mean_qual = sum(record.letter_annotations.get("phred_quality", [])) / seq_len if seq_len else 0

            if (gc_min <= seq_gc <= gc_max and
                len_min <= seq_len <= len_max and
                mean_qual >= quality_threshold):

                passed_count += 1
                filtered[record.id] = record
                SeqIO.write(record, out_handle, "fastq")

    logging.info(f"Finished filtering. Total: {total_count}, Passed: {passed_count}")

    return {
        "total": total_count,
        "passed": passed_count,
        "filtered_records": filtered
    }

# ======================================
# CLI
# ======================================
def parse_args():
    parser = argparse.ArgumentParser(description="FASTQ filtering tool")

    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output FASTQ file")

    parser.add_argument("--gc-min", type=float, default=0)
    parser.add_argument("--gc-max", type=float, default=100)

    parser.add_argument("--len-min", type=int, default=0)
    parser.add_argument("--len-max", type=int, default=int(1e9))

    parser.add_argument("--qual", type=float, default=0)

    return parser.parse_args()


def main():
    args = parse_args()

    result = filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=(args.gc_min, args.gc_max),
        length_bounds=(args.len_min, args.len_max),
        quality_threshold=args.qual
    )

    print(f"Total reads: {result['total']}")
    print(f"Passed reads: {result['passed']}")


if __name__ == "__main__":
    main()