"""Consensus, conservation, and variant-column analysis of an alignment."""

from collections import Counter


def _columns(aligned):
    return zip(*aligned) if aligned else []


def consensus(aligned):
    """Most common residue per column (deterministic tie-break by residue)."""
    out = []
    for col in _columns(aligned):
        counts = Counter(col)
        top = max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0]
        out.append(top)
    return "".join(out)


def conservation(aligned):
    """Fraction of rows matching the column consensus, per column."""
    cons = consensus(aligned)
    rows = len(aligned)
    return [sum(1 for c in col if c == cons[i]) / rows for i, col in enumerate(_columns(aligned))]


def variant_columns(aligned):
    """Indices of columns where not all rows agree."""
    return [i for i, col in enumerate(_columns(aligned)) if len(set(col)) > 1]
