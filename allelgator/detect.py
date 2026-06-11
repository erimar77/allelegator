"""Classify a sequence as protein or nucleotide."""

_NUC = set("ACGTUN")


def detect_type(seq):
    """Return 'nucleotide' if the sequence is overwhelmingly ACGTUN, else 'protein'."""
    letters = [c for c in seq.upper() if c.isalpha()]
    if not letters:
        return "protein"
    nuc = sum(1 for c in letters if c in _NUC)
    return "nucleotide" if nuc / len(letters) > 0.85 else "protein"
