"""Render alignments to the terminal with rich: full (T-COFFEE-style) and compact diff views."""

from rich.console import Group
from rich.table import Table
from rich.text import Text

_MAX_LABEL = 15


def _labeler(names):
    """Return (width, fn) where fn(name) -> padded, truncated label."""
    longest = max([len(n) for n in names] + [len("cons")])
    width = min(longest, _MAX_LABEL) + 1

    def fn(name):
        return name[: width - 1].ljust(width)

    return width, fn


def render_full(names, aligned, variant_cols, start=1, width=None,
                annotations=None, show_consensus=True):
    """Build a rich renderable for the full alignment view.

    Wrapped into blocks with a numbered ruler. Variant columns are bold red
    across every row; a `cons` row marks fully-conserved columns with `*`
    (T-COFFEE style). `annotations` is reserved for a future boxed feature
    track and is currently unused.
    """
    variant = set(variant_cols)
    total = len(aligned[0]) if aligned else 0
    block = width or 60
    lbl, label = _labeler(names)
    cons = ("".join("*" if len(set(col)) == 1 else " " for col in zip(*aligned))
            if aligned else "")

    lines = []
    for off in range(0, total, block):
        cols = list(range(off, min(off + block, total)))
        # Ruler: start residue number on the left, end number on the right.
        left = str(start + off)
        right = str(start + cols[-1])
        ruler = Text(" " * lbl + left, style="bold cyan")
        pad = (lbl + len(cols)) - len(ruler.plain) - len(right)
        ruler.append(" " * max(1, pad) + right, style="bold cyan")
        lines.append(ruler)
        # One row per allele.
        for name, seq in zip(names, aligned):
            row = Text(label(name), style="bold")
            for c in cols:
                ch = seq[c]
                if c in variant:
                    style = "bold red"
                elif ch == "-":
                    style = "dim"
                else:
                    style = ""
                row.append(ch, style=style)
            lines.append(row)
        # Consensus row.
        if show_consensus:
            crow = Text(label("cons"), style="dim")
            for c in cols:
                crow.append(cons[c], style="green" if cons[c] == "*" else "")
            lines.append(crow)
        lines.append(Text(""))
    return Group(*lines)


def render_reference(names, aligned, variant_cols, ref_index=0, start=1,
                     width=None, show_consensus=True):
    """Dot-identity view: one reference row shown in full; every other row shows
    `.` where it matches the reference and the differing residue otherwise.

    The consensus row marks fully-conserved columns with `*` and the rest with
    `.` (the convention used by published dot-identity alignments).
    """
    total = len(aligned[0]) if aligned else 0
    block = width or 60
    lbl, label = _labeler(names)
    ref = aligned[ref_index] if aligned else ""
    cons = ("".join("*" if len(set(col)) == 1 else "." for col in zip(*aligned))
            if aligned else "")

    lines = []
    for off in range(0, total, block):
        cols = list(range(off, min(off + block, total)))
        left = str(start + off)
        right = str(start + cols[-1])
        ruler = Text(" " * lbl + left, style="bold cyan")
        pad = (lbl + len(cols)) - len(ruler.plain) - len(right)
        ruler.append(" " * max(1, pad) + right, style="bold cyan")
        lines.append(ruler)
        for idx, (name, seq) in enumerate(zip(names, aligned)):
            row = Text(label(name), style="bold")
            for c in cols:
                ch = seq[c]
                if idx == ref_index:
                    row.append(ch, style="dim" if ch == "-" else "")
                elif ch == ref[c]:
                    row.append(".", style="dim")
                else:
                    row.append(ch, style="bold red")
            lines.append(row)
        if show_consensus:
            crow = Text(label("cons"), style="dim")
            for c in cols:
                crow.append(cons[c], style="green" if cons[c] == "*" else "dim")
            lines.append(crow)
        lines.append(Text(""))
    return Group(*lines)


def render_diff(names, aligned, variant_cols, start=1):
    """Compact view: a summary line plus one table row per variant column."""
    total = len(aligned[0]) if aligned else 0
    pct = 100.0 * (total - len(variant_cols)) / total if total else 100.0
    summary = Text(
        f"Alignment length: {total}   "
        f"Identical columns: {pct:.1f}%   "
        f"Variant columns: {len(variant_cols)}",
        style="bold",
    )
    table = Table(show_edge=True, header_style="bold cyan")
    table.add_column("pos", justify="right")
    for n in names:
        table.add_column(n[:_MAX_LABEL], justify="center")
    for c in variant_cols:
        cells = [str(start + c)]
        cells.extend(seq[c] for seq in aligned)
        table.add_row(*cells)
    return Group(summary, Text(""), table)
