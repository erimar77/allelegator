"""Command-line entry point: gather sequences, align, and render."""

import argparse
import sys
from pathlib import Path

from rich.console import Console

from . import __version__, align, analyze, detect, loaders, render


def _split_label(token):
    """Split an `ID=Label` token into (source, label); label is None if absent."""
    if "=" in token:
        source, label = token.split("=", 1)
        return source, label
    return token, None


def _relabel(records, label):
    """Apply a friendly label to loaded records (suffixing if there are several)."""
    if label is None:
        return records
    if len(records) == 1:
        return [(label, records[0][1])]
    return [(f"{label}{i + 1}", seq) for i, (_, seq) in enumerate(records)]


def _gather(args):
    """Return a list of (name, seq). Raises loaders.LoaderError on failure.

    Any input token may carry a friendly name as `SOURCE=Label`, e.g.
    `NP_416871.1=MG1655` or `alleles.fasta=MyGene`.
    """
    if args.sample:
        return loaders.load_sample()
    records = []
    for token in args.inputs + (args.fetch or []):
        source, label = _split_label(token)
        if Path(source).exists():
            recs = loaders.parse_fasta(Path(source).read_text())
        else:
            recs = [loaders.fetch_ncbi(source, db=args.db)]
        records.extend(_relabel(recs, label))
    return records


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="allelegator",
        description="AlleleGator 🐊 — line up your alleles and see what bites.",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show version and exit",
    )
    p.add_argument("inputs", nargs="*", help="FASTA files and/or NCBI accessions")
    p.add_argument("--sample", action="store_true", help="use the bundled APOE sample")
    p.add_argument("--diff", action="store_true", help="show compact variant-only view")
    p.add_argument(
        "--ref",
        nargs="?",
        const="",
        default=None,
        metavar="NAME",
        help="dot-identity view against a reference (default: first sequence)",
    )
    p.add_argument("--fetch", nargs="+", help="explicitly fetch these NCBI accessions")
    p.add_argument("--db", help="force NCBI db (protein|nuccore)")
    p.add_argument("--start", type=int, default=1, help="first residue number for the ruler")
    p.add_argument("--width", type=int, default=None, help="columns per block")
    p.add_argument("--no-color", action="store_true", help="disable color")
    args = p.parse_args(argv)

    console = Console(no_color=args.no_color)

    if not args.sample and not args.inputs and not args.fetch:
        console.print("[red]error:[/] no inputs. Pass FASTA files/accessions or --sample.")
        return 2

    try:
        records = _gather(args)
    except loaders.LoaderError as e:
        console.print(f"[red]error:[/] {e}")
        return 1

    if len(records) < 2:
        console.print("[yellow]warning:[/] need at least two sequences to align.")
        return 1

    names = [n for n, _ in records]
    seqs = [s for _, s in records]
    types = {detect.detect_type(s) for s in seqs}
    if len(types) > 1:
        console.print("[red]error:[/] mixed protein and nucleotide inputs.")
        return 1
    seq_type = types.pop()

    aligned = align.center_star_msa(seqs, seq_type)
    variant = analyze.variant_columns(aligned)

    if args.diff:
        console.print(render.render_diff(names, aligned, variant, start=args.start))
    elif args.ref is not None:
        ref_name = args.ref or names[0]
        if ref_name not in names:
            console.print(
                f"[red]error:[/] reference '{ref_name}' not found. "
                f"Available: {', '.join(names)}"
            )
            return 1
        console.print(render.render_reference(
            names, aligned, variant, ref_index=names.index(ref_name),
            start=args.start, width=args.width,
        ))
    else:
        console.print(render.render_full(names, aligned, variant, start=args.start, width=args.width))
    return 0


if __name__ == "__main__":
    sys.exit(main())
