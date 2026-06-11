"""Load allele sequences from FASTA text, local files, NCBI, or the bundled sample."""

import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

_SAMPLE = Path(__file__).parent / "data" / "apoe_alleles.fasta"
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class LoaderError(Exception):
    """Raised when sequences cannot be loaded."""


def parse_fasta(text):
    """Parse FASTA text into a list of (name, sequence) tuples."""
    records = []
    name = None
    chunks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            parts = line[1:].split()
            name = parts[0] if parts else ""
            chunks = []
        else:
            chunks.append("".join(line.split()).upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise LoaderError("no FASTA records found")
    return records


def load_sample():
    """Load the bundled APOE allele sample (ε3/ε2/ε4)."""
    return parse_fasta(_SAMPLE.read_text())


def fetch_ncbi(accession, db=None):
    """Fetch a single FASTA record from NCBI E-utilities by accession.

    Returns (accession, sequence). If `db` is not given, the database is guessed
    from the accession and the other database is tried as a fallback. Raises
    LoaderError on network/parse failure.
    """
    candidates = [db] if db else _candidate_dbs(accession)
    last_err = None
    for d in candidates:
        try:
            text = _efetch(accession, d)
        except OSError as e:
            last_err = e
            continue
        records = parse_fasta_or_empty(text)
        if records:
            return (accession, records[0][1])
        last_err = last_err or f"no sequence in db={d}"
    raise LoaderError(f"NCBI fetch failed for {accession}: {last_err}")


def _efetch(accession, db):
    params = urlencode({"db": db, "id": accession, "rettype": "fasta", "retmode": "text"})
    with urlopen(f"{_EFETCH}?{params}", timeout=30) as resp:
        return resp.read().decode()


def parse_fasta_or_empty(text):
    """Like parse_fasta but returns [] instead of raising on no records."""
    try:
        return parse_fasta(text)
    except LoaderError:
        return []


def _candidate_dbs(accession):
    """Guess the most likely NCBI database, with the other as a fallback."""
    primary = _guess_db(accession)
    return [primary, "nuccore" if primary == "protein" else "protein"]


def _guess_db(accession):
    """Heuristic protein-vs-nucleotide guess from an accession.

    RefSeq protein prefixes (NP_/XP_/YP_/WP_/AP_) and GenBank protein accessions
    (three leading letters, e.g. AAN81356, WPV06440) are protein; one- or
    two-letter GenBank prefixes and other RefSeq prefixes are nucleotide.
    """
    acc = accession.split(".")[0]
    if "_" in acc:
        return "protein" if acc[:2] in ("NP", "XP", "YP", "WP", "AP") else "nuccore"
    m = re.match(r"^([A-Za-z]+)", acc)
    return "protein" if (m and len(m.group(1)) >= 3) else "nuccore"
