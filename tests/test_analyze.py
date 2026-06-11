from allelgator.analyze import conservation, consensus, variant_columns

ALN = ["ACGT", "ACGA", "ACGT"]


def test_consensus():
    assert consensus(ALN) == "ACGT"


def test_conservation():
    assert conservation(ALN) == [1.0, 1.0, 1.0, 2 / 3]


def test_variant_columns():
    assert variant_columns(ALN) == [3]


def test_no_variants():
    assert variant_columns(["AAA", "AAA"]) == []
