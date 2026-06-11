from rich.console import Console

from allelgator.render import render_diff, render_full


def _plain(renderable):
    con = Console(width=120, no_color=True)
    with con.capture() as cap:
        con.print(renderable)
    return cap.get()


def test_full_view_contains_labels_and_ruler():
    names = ["e3", "e2", "e4"]
    aligned = ["ACGTACGT", "ACGTTCGT", "ACGTACGT"]
    out = _plain(render_full(names, aligned, variant_cols=[4], start=1))
    assert "e3" in out and "e2" in out and "e4" in out
    assert "1" in out


def test_full_view_has_consensus_row():
    names = ["a", "b"]
    aligned = ["ACGT", "ACGA"]
    out = _plain(render_full(names, aligned, variant_cols=[3], start=1))
    assert "cons" in out
    assert "*" in out  # conserved columns marked


def test_full_view_wraps_blocks():
    names = ["a", "b"]
    aligned = ["A" * 200, "A" * 200]
    out = _plain(render_full(names, aligned, variant_cols=[], start=1, width=50))
    assert out.count("\n") > 4


def test_diff_view_lists_variant_positions():
    names = ["e3", "e2", "e4"]
    aligned = ["ACGTACGT", "ACGTTCGT", "ACGAACGT"]
    out = _plain(render_diff(names, aligned, variant_cols=[3, 4], start=1))
    assert "4" in out and "5" in out
    assert "%" in out
