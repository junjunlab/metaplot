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


def CalculatePeaksRelativeDist(peaksfa, motif=DEFAULT_MOTIF, mid=3):
    """Distance of every motif hit to the centre of its peak sequence.

    ``mid`` is added to the 0-based start of each match to pick the base of
    interest inside the motif; the default (3) targets the A of "RRACH" when
    the motif is matched in its reverse-complement orientation.
    """
    pattern = re.compile(motif, flags=re.IGNORECASE)
    offset = int(mid)

    relpos = []
    for _, seq in iter_fasta(peaksfa):
        if not seq:
            continue
        # unchanged from the original implementation: integer-ish centre
        peakmid = round(len(seq) / 2, ndigits=0)
        for match in pattern.finditer(seq):
            relpos.append(match.start() + offset - peakmid)

    return relpos
