# AlleleGator 🐊

**Line up your alleles and see what bites.**

A terminal multiple-sequence aligner and visualizer. Point it at several allele
sequences of one gene — from local FASTA files, NCBI accessions, or the bundled
sample — and it renders a colored alignment where the columns that *differ* jump
out. The whole point: a handful of allele sequences is often only ~100 bytes of
text, yet impossible to compare by eye in its raw form. AlleleGator fixes that.

Pure Python, one dependency (`rich`). No Biopython, no MUSCLE/MAFFT, no build step.

---

## Install

```bash
python -m venv .venv
.venv/bin/pip install -e .
```

This puts the command `allelegator` on your PATH (inside the venv), plus a
shorter alias `gator` — the two are identical.

## Quick start

```bash
# Bundled demo: the three APOE alleles (ε2/ε3/ε4)
allelegator --sample

# The "just tell me what differs" view
allelegator --sample --diff

# Your own sequences
allelegator alleles.fasta
allelegator seq_a.fasta seq_b.fasta seq_c.fasta

# Pull straight from NCBI by accession (the way the data is shared)
allelegator NP_000032.1 NP_000041.1

# Strain comparison: EvgS sensor kinase across E. coli strains, dot-identity
# view against the K-12 (MG1655) reference, fetched live from NCBI.
# Use SOURCE=Label to give each accession a friendly name.
allelegator NP_416871.1=MG1655 WPV06440.1=O157 AAN81356.1=CFT073 \
      --ref MG1655
```

Any input — a file path or an accession — may carry a friendly name as
`SOURCE=Label` (e.g. `NP_416871.1=MG1655`). The label is what shows up in the
alignment and what you pass to `--ref`.

The protein vs. nucleotide database is guessed from each accession (with the
other database tried as a fallback), so a mix of RefSeq and GenBank protein IDs
just works.

## What the output looks like

**Full view** (T-COFFEE-style): wrapped blocks with a numbered ruler, one row
per allele, and a `cons` row marking fully-conserved columns with `*`. Columns
where the alleles disagree are shown in **bold red** across every row.

```
     61                                                       120
e3   RALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVCGRLVQYRG
e2   RALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVCGRLVQYRG
e4   RALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVRGRLVQYRG
cons *************************************************** ********
```

**Reference / dot-identity view** (`--ref [NAME]`): one reference sequence (the
top row, or the one you name) is shown in full; every other row shows a `.`
where it matches the reference and the actual residue only where it differs.
This is the "inverted" style common in strain-comparison figures — blanks mean
"same," letters mean "changed."

```
            1                                                         60
NP_416871.1 MKFLPYIFLLCCGLWSTISFADEDYIEYRGISSNNRVTLDPLRLSNKELRWLASKKNLVI
WPV06440.1  ......................G.....................................
AAN81356.1  ...........X..........G.....H...............................
CBJ01993.1  ............A....M..........................................
cons        ***********..****.****.*****.*******************************
```

**Diff view** (`--diff`): a one-line summary plus a table of just the variant
columns — for highly similar alleles this is usually all you want.

```
Alignment length: 299   Identical columns: 99.3%   Variant columns: 2

  pos    e3   e2   e4
  112    C    C    R
  158    R    C    R
```

## Usage

```
allelegator [INPUTS...] [options]

INPUTS                Zero or more FASTA file paths and/or NCBI accession IDs
                      (mixed is fine). Anything that exists as a file is read as
                      FASTA; everything else is treated as an accession to fetch.

  --sample            Use the bundled APOE allele sample (ignores INPUTS).
  --diff              Show the compact variant-only view.
  --ref [NAME]        Dot-identity view against a reference (default: first
                      sequence). Matches show as `.`, differences as the residue.
  --fetch ID [ID...]  Explicitly fetch these accessions from NCBI.
  --db DB             Force the NCBI database (protein | nuccore).
  --start N           First residue/base number shown on the ruler (default 1).
  --width N           Columns per wrapped block (default 60).
  --no-color          Disable color (for piping or plain terminals).
  -h, --help          Show help.
```

Protein vs. nucleotide input is auto-detected; mixing the two in one run is an
error.

## How it aligns

AlleleGator builds a **center-star multiple sequence alignment**: it aligns
every sequence to the most central one (highest summed pairwise similarity)
using Needleman–Wunsch global alignment, then merges induced gaps so columns
line up. The gap penalty is tuned harsher than a mismatch, so near-identical
alleles align by substitution instead of sprouting spurious gaps. This targets
the "alleles of one gene" case — it is not meant to replace MUSCLE/MAFFT on
highly divergent sequences.

Validated against a T-COFFEE alignment of *E. coli* `gadA`/`gadB` (1401 columns,
31 variant positions, zero gaps) — the conserved-column pattern matches exactly.
See `examples/ecoli_gad.fasta`.

## Tests

```bash
.venv/bin/pytest
```

## Man page

```bash
man -l man/allelegator.1
```

## License

MIT.
