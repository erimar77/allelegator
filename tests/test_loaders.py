import pytest

from allelgator import loaders
from allelgator.align import center_star_msa
from allelgator.analyze import variant_columns
from allelgator.loaders import LoaderError, fetch_ncbi, load_sample, parse_fasta


def test_parse_two_records():
    text = ">e3 apoe\nHLRK\nLRKR\n>e4\nHLRKLRKR\n"
    assert parse_fasta(text) == [("e3", "HLRKLRKR"), ("e4", "HLRKLRKR")]


def test_parse_ignores_blank_lines_and_uppercases():
    assert parse_fasta(">a\nhl rk\n\n") == [("a", "HLRK")]


def test_parse_keeps_piped_names():
    assert parse_fasta(">gnl|ECOLI|EG114 gadB\nACGT\n") == [("gnl|ECOLI|EG114", "ACGT")]


def test_parse_empty_raises():
    with pytest.raises(LoaderError):
        parse_fasta("   \n")


def test_sample_has_three_alleles():
    assert [n for n, _ in load_sample()] == ["e3", "e2", "e4"]


def test_sample_has_exactly_two_variant_columns():
    recs = load_sample()
    aligned = center_star_msa([s for _, s in recs], "protein")
    cols = variant_columns(aligned)
    assert len(cols) == 2, f"expected 2 variant columns, got {cols}"


def test_fetch_ncbi_parses_fasta(monkeypatch):
    fake = ">NP_000032.1 apoe\nHLRKLRKR\n"

    class FakeResp:
        def read(self):
            return fake.encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(loaders, "urlopen", lambda url, timeout=0: FakeResp())
    name, seq = fetch_ncbi("NP_000032.1", db="protein")
    assert name == "NP_000032.1"
    assert seq == "HLRKLRKR"


def test_guess_db():
    from allelgator.loaders import _guess_db

    assert _guess_db("NP_416871.1") == "protein"
    assert _guess_db("WPV06440.1") == "protein"   # 3-letter GenBank protein
    assert _guess_db("AAN81356.1") == "protein"
    assert _guess_db("NM_000041.3") == "nuccore"
    assert _guess_db("U00096.3") == "nuccore"      # 1-letter GenBank nucleotide


def test_fetch_ncbi_falls_back_to_other_db(monkeypatch):
    # First db (guessed) 400s; fallback db succeeds.
    calls = []

    def fake_efetch(accession, db):
        calls.append(db)
        if db == calls[0]:
            raise OSError("HTTP Error 400")
        return ">x\nMKLV\n"

    monkeypatch.setattr(loaders, "_efetch", fake_efetch)
    name, seq = fetch_ncbi("X12345.1")
    assert seq == "MKLV"
    assert len(calls) == 2  # tried primary then fallback


def test_fetch_ncbi_network_error(monkeypatch):
    def boom(url, timeout=0):
        raise OSError("no network")

    monkeypatch.setattr(loaders, "urlopen", boom)
    with pytest.raises(LoaderError):
        fetch_ncbi("NP_000032.1", db="protein")
