"""Command-line entry point: gather sequences, align, and render."""

import argparse
import sys
from pathlib import Path

from rich.console import Console

from . import align, analyze, detect, loaders, render


def _gather(args):
    """Return a list of (name, seq). Raises loaders.LoaderError on failure."""
    if args.sample:
        return loaders.load_sample()
    records = []
    for token in args.inputs + (args.fetch or []):
        if Path(token).exists():
            records.extend(loaders.parse_fasta(Path(token).read_text()))
        else:
            records.append(loaders.fetch_ncbi(token, db=args.db))
    return records


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="gator",
        description="AlleleGator 🐊 — line up your alleles and see what bites.",
    )
    p.add_argument("inputs", nargs="*", help="FASTA files and/or NCBI accessions")
    p.add_argument("--sample", action="store_true", help="use the bundled APOE sample")
    p.add_argument("--diff", action="store_true", help="show compact variant-only view")
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
    else:
        console.print(render.render_full(names, aligned, variant, start=args.start, width=args.width))
    return 0


if __name__ == "__main__":
    sys.exit(main())
