from allelgator.cli import main


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
