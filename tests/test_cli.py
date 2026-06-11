from allelgator.cli import _relabel, _split_label, main


def test_split_label():
    assert _split_label("NP_416871.1=MG1655") == ("NP_416871.1", "MG1655")
    assert _split_label("alleles.fasta") == ("alleles.fasta", None)


def test_relabel_single_record():
    assert _relabel([("acc", "MKLV")], "MG1655") == [("MG1655", "MKLV")]


def test_relabel_multi_record_suffixes():
    recs = [("a", "AA"), ("b", "CC")]
    assert _relabel(recs, "Gene") == [("Gene1", "AA"), ("Gene2", "CC")]


def test_relabel_none_passthrough():
    recs = [("a", "AA")]
    assert _relabel(recs, None) == recs


def test_sample_full_view(capsys):
    rc = main(["--sample"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "e3" in out and "e2" in out and "e4" in out


def test_sample_diff_view(capsys):
    rc = main(["--sample", "--diff"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "Variant columns: 2" in out


def test_no_input_errors(capsys):
    rc = main([])
    assert rc != 0
