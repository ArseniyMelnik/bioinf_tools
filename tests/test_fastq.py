import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
import os
from main import DNASequence, RNASequence, AminoAcidSequence, filter_fastq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def test_dna_complement():
    dna = DNASequence("ATGC")
    assert str(dna.complement()) == "DNASequence: TACG"


def test_rna_reverse_complement():
    rna = RNASequence("AUGC")
    assert str(rna.reverse_complement()) == "RNASequence: GCAU"


def test_invalid_dna():
    with pytest.raises(ValueError):
        DNASequence("ATBX")


def test_aa_composition():
    protein = AminoAcidSequence("AAA")
    comp = protein.aa_composition()
    assert comp["A"] == 100.0


def test_filter_gc(tmp_path):
    input_file = tmp_path / "test.fastq"
    output_file = tmp_path / "out.fastq"

    records = [
        SeqRecord(Seq("GGGG"), id="high_gc", letter_annotations={"phred_quality":[40]*4}),
        SeqRecord(Seq("ATAT"), id="low_gc", letter_annotations={"phred_quality":[40]*4})
    ]

    SeqIO.write(records, input_file, "fastq")

    result = filter_fastq(str(input_file), str(output_file), gc_bounds=(50,100))

    assert result["passed"] == 1


def test_filter_length(tmp_path):
    input_file = tmp_path / "test.fastq"
    output_file = tmp_path / "out.fastq"

    records = [
        SeqRecord(Seq("AT"), id="short", letter_annotations={"phred_quality":[40]*2}),
        SeqRecord(Seq("ATGCGA"), id="long", letter_annotations={"phred_quality":[40]*6})
    ]

    SeqIO.write(records, input_file, "fastq")

    result = filter_fastq(str(input_file), str(output_file), length_bounds=(5,10))

    assert result["passed"] == 1


def test_filter_quality(tmp_path):
    input_file = tmp_path / "test.fastq"
    output_file = tmp_path / "out.fastq"

    records = [
        SeqRecord(Seq("ATGC"), id="low_q", letter_annotations={"phred_quality":[10]*4}),
        SeqRecord(Seq("ATGC"), id="high_q", letter_annotations={"phred_quality":[40]*4})
    ]

    SeqIO.write(records, input_file, "fastq")

    result = filter_fastq(str(input_file), str(output_file), quality_threshold=30)

    assert result["passed"] == 1


def test_output_file_created(tmp_path):
    input_file = tmp_path / "test.fastq"
    output_file = tmp_path / "out.fastq"

    records = [
        SeqRecord(Seq("ATGC"), id="r1", letter_annotations={"phred_quality":[40]*4})
    ]

    SeqIO.write(records, input_file, "fastq")

    filter_fastq(str(input_file), str(output_file))

    assert os.path.exists(output_file)
