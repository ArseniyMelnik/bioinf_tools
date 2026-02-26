# bio_files_processor.py

import os

# -------------------------------
# 1) Convert multiline FASTA to single-line FASTA
# -------------------------------
def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Reads a multi-line FASTA file and writes a new FASTA file
    where each sequence is in a single line.

    Parameters:
        input_fasta: path to the input FASTA file
        output_fasta: path to output FASTA file (optional)
                      if None, saves as 'oneline_<input_fasta>'
    """
    if output_fasta is None:
        # Construct default output filename
        base_name = os.path.basename(input_fasta)
        output_fasta = os.path.join(
            os.path.dirname(input_fasta), f'oneline_{base_name}'
        )

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        seq_name = ''
        seq_lines = []
        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                # If there is a previous sequence, write it to output
                if seq_name:
                    outfile.write(seq_name + ''.join(seq_lines) + '\n')
                seq_name = line + '\n'
                seq_lines = []
            else:
                seq_lines.append(line)
        # Write the last sequence
        if seq_name:
            outfile.write(seq_name + ''.join(seq_lines) + '\n')


# -------------------------------
# 2) Parse BLAST output
# -------------------------------
def parse_blast_output(input_file: str, output_file: str):
    """
    Parses a BLAST output file (txt format) to extract the top hit
    description for each query and writes a sorted list of hits.

    Parameters:
        input_file: path to the BLAST txt output
        output_file: path to save extracted protein names
    """
    hits = set()
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        # Look for query blocks
        if line.startswith('Sequences producing significant alignments:'):
            i += 1  # Move to first hit
            while i < len(lines) and lines[i].strip() == '':
                i += 1
            if i < len(lines):
                # Take the first non-empty hit line
                hit_line = lines[i].strip()
                # The Description is assumed to be the full line after whitespace
                if hit_line:
                    # split and take all after first whitespace (simplest)
                    parts = hit_line.split(None, 1)
                    if len(parts) == 2:
                        hits.add(parts[1])
            i += 1
        else:
            i += 1

    # Write sorted unique hits to output
    with open(output_file, 'w') as outfile:
        for h in sorted(hits):
            outfile.write(h + '\n')
