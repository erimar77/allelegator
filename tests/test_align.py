from allelgator.align import center_star_msa, needleman_wunsch, nucleotide_score


def test_identical_sequences_no_gaps():
    a, b, score = needleman_wunsch("ACGT", "ACGT", nucleotide_score, gap=-1)
    assert a == "ACGT"
    assert b == "ACGT"
    assert score == 4


def test_single_insertion():
    a, b, score = needleman_wunsch("ACGT", "ACGGT", nucleotide_score, gap=-1)
    assert len(a) == len(b)
    assert "-" in a
    assert a.replace("-", "") == "ACGT"
    assert b.replace("-", "") == "ACGGT"


def test_msa_equal_length_and_order():
    seqs = ["ACGTACGT", "ACGTTACGT", "ACGTACGT"]
    out = center_star_msa(seqs, "nucleotide")
    assert len({len(s) for s in out}) == 1
    assert len(out) == 3
    assert [s.replace("-", "") for s in out] == seqs


def test_msa_single_sequence():
    assert center_star_msa(["ACGT"], "nucleotide") == ["ACGT"]


def test_msa_identical_no_gaps():
    assert center_star_msa(["AAA", "AAA"], "nucleotide") == ["AAA", "AAA"]


def test_msa_substitutions_only_no_gaps():
    # Equal-length sequences differing only by substitutions must not gain gaps
    # (regression: gap penalty must beat a mismatch).
    a = "ATGGATAAGAAGCAAGTAACGGATTTAAGG"
    b = "ATGGACCAGAAGCTGTTAACGGATTTCCGC"
    out = center_star_msa([a, b], "nucleotide")
    assert all("-" not in row for row in out)
    assert out == [a, b]
