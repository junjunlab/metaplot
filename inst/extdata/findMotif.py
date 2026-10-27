import re

# default motif: "RRACH" plus its reverse complement "DGGYY"
DEFAULT_MOTIF = r"[GA][GA]AC[ACT]|[TGA]GT[CT][CT]"


def iter_fasta(path):
    """Yield ``(header, sequence)`` pairs from a FASTA file.

    Sequences wrapped over several lines are joined back into one string;
    a header line is never treated as sequence data.
    """
    if path is None:
        raise ValueError("a FASTA file path must be supplied")

    name = None
    chunks = []

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:]
                chunks = []
            elif name is None:
                # sequence before the first header: keep it under a dummy name
                name = ""
                chunks.append(line)
            else:
                chunks.append(line)

    if name is not None:
        yield name, "".join(chunks)


def findMotif(peaksfa, motif=DEFAULT_MOTIF):
    """Return every motif occurrence found in the FASTA file."""
    pattern = re.compile(motif, flags=re.IGNORECASE)

    hits = []
    for _, seq in iter_fasta(peaksfa):
        hits.extend(pattern.findall(seq))

    return hits
