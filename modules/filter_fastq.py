# modules/filter_fastq.py

import os

def calc_gc(seq: str) -> float:
    """Calculate GC content as percentage."""
    return 100 * sum(1 for b in seq.upper() if b in "GC") / len(seq)

def calc_mean_quality(qual: str) -> float:
    """Convert Phred33 quality string to numeric scores and return mean quality."""
    if not qual:
        return 0
    scores = [ord(ch) - 33 for ch in qual]
    return sum(scores) / len(scores)

def check_bounds(value, bounds) -> bool:
    """Check if value is within bounds."""
    if isinstance(bounds, (int, float)):
        lower, upper = 0, bounds
    else:
        lower, upper = bounds
    return lower <= value <= upper

def filter_fastq(input_fastq: str,
                 output_fastq: str,
                 gc_bounds=(0,100),
                 length_bounds=(0,2**32),
                 quality_threshold=0) -> dict:
    """
    Filter a FASTQ file on-the-fly.
    Reads FASTQ file, applies filters, writes passing sequences to output,
    and returns a dictionary with filtered sequences.

    Parameters:
        input_fastq: path to input FASTQ
        output_fastq: path to output FASTQ (will be created if missing)
        gc_bounds: tuple or int
        length_bounds: tuple or int
        quality_threshold: minimum mean Phred33 quality

    Returns:
        dict: {sequence_name: (sequence, quality_string)}
    """
    # Ensure output directory exists
    out_dir = os.path.dirname(output_fastq)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    filtered = {}

    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        while True:
            # Read 4 lines: name, seq, +, quality
            name_line = infile.readline()
            if not name_line:
                break  # EOF
            seq_line = infile.readline().rstrip()
            plus_line = infile.readline()
            qual_line = infile.readline().rstrip()

            # Apply filters
            if (check_bounds(calc_gc(seq_line), gc_bounds) and
                check_bounds(len(seq_line), length_bounds) and
                calc_mean_quality(qual_line) >= quality_threshold):

                # Add to dictionary
                filtered[name_line.strip()] = (seq_line, qual_line)

                # Write to output file
                outfile.write(name_line)
                outfile.write(seq_line + '\n')
                outfile.write(plus_line)
                outfile.write(qual_line + '\n')

    return filtered