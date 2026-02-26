# modules/filter_fastq.py


def calc_gc(seq):
    """Calculate GC content as percentage."""
    if not seq:
        return 0
    return 100 * sum(1 for b in seq.upper() if b in "GC") / len(seq)


def calc_length(seq):
    """Return sequence length."""
    return len(seq)


def phred33_to_q(qual):
    """Convert Phred33 quality string to numeric scores."""
    return [ord(ch) - 33 for ch in qual]


def calc_mean_quality(qual):
    """Return mean quality score."""
    scores = phred33_to_q(qual)
    return sum(scores) / len(scores) if scores else 0


def check_bounds(value, bounds):
    """
    Check if value is within bounds.
    If bounds is a single number, treat it as upper bound, lower = 0.
    """
    if isinstance(bounds, (int, float)):
        lower, upper = 0, bounds
    else:
        lower, upper = bounds
    return lower <= value <= upper


def filter_fastq(
    seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0
):
    """
    Filter a dictionary of FASTQ sequences by GC, length, and mean quality.

    Parameters:
        seqs: dict {name: (sequence, quality_string)}
        gc_bounds: tuple(min, max) or single number (upper bound)
        length_bounds: tuple(min, max) or single number (upper bound)
        quality_threshold: minimum average Phred33 quality

    Returns:
        dict with filtered sequences
    """
    filtered = {}
    for name, (seq, qual) in seqs.items():
        if (
            check_bounds(calc_gc(seq), gc_bounds)
            and check_bounds(calc_length(seq), length_bounds)
            and calc_mean_quality(qual) >= quality_threshold
        ):
            filtered[name] = (seq, qual)
    return filtered
