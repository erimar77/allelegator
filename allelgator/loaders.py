"""Load allele sequences from FASTA text, local files, NCBI, or the bundled sample."""

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

    Returns (accession, sequence). Raises LoaderError on network/parse failure.
    """
    db = db or _guess_db(accession)
    params = urlencode({"db": db, "id": accession, "rettype": "fasta", "retmode": "text"})
    url = f"{_EFETCH}?{params}"
    try:
        with urlopen(url, timeout=30) as resp:
            text = resp.read().decode()
    except OSError as e:
        raise LoaderError(f"NCBI fetch failed for {accession}: {e}") from e
    records = parse_fasta(text)
    if not records:
        raise LoaderError(f"NCBI returned no sequence for {accession}")
    return (accession, records[0][1])


def _guess_db(accession):
    """Protein RefSeq accessions start with NP_/XP_/AP_/YP_; default to nuccore."""
    return "protein" if accession[:3] in ("NP_", "XP_", "AP_", "YP_") else "nuccore"
