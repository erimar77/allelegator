from allelgator.detect import detect_type


def test_detects_nucleotide():
    assert detect_type("ACGTACGTACGTN") == "nucleotide"


def test_detects_protein():
    assert detect_type("HLRKLRKRLLRDADDLQKR") == "protein"


def test_gaps_ignored():
    assert detect_type("ACGT--ACGT") == "nucleotide"
