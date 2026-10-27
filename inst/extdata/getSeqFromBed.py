from pyfaidx import Fasta

VALID_TYPES = ("bed", "narrowpeak")


def GetPeaksFa(peak, genomefie, outfasta, type='bed', a=0, b=0, minlen=5):
    """Write the sequence of every peak in ``peak`` to ``outfasta``.

    ``a`` and ``b`` extend the peak upstream and downstream **in transcript
    orientation**: on the minus strand "upstream" means towards higher
    coordinates. Both BED and narrowPeak are 0-based half-open, so no
    coordinate adjustment is needed for either format.

    Returns the number of sequences written.
    """
    if peak is None or genomefie is None or outfasta is None:
        raise ValueError("`peak`, `genomefie` and `outfasta` must all be supplied.")
    if type not in VALID_TYPES:
        raise ValueError(
            "`type` must be one of %s, got %r" % (", ".join(VALID_TYPES), type)
        )

    a = int(a)
    b = int(b)
    minlen = int(minlen)

    genome = Fasta(genomefie)
    written = 0

    with open(peak, 'r') as bed, open(outfasta, 'w') as outfa:
        for line in bed:
            line = line.strip()
            if not line or line.startswith(('#', 'track', 'browser')):
                continue

            fields = line.split()
            if len(fields) < 6:
                raise ValueError(
                    "each peak needs at least 6 BED columns, got %d: %s"
                    % (len(fields), line)
                )

            chrom = fields[0]
            strand = fields[5]
            start = int(fields[1])
            end = int(fields[2])

            if chrom not in genome:
                raise KeyError("sequence %r not found in %r" % (chrom, genomefie))

            # BED and narrowPeak are both 0-based half-open
            if strand == '-':
                # upstream (a) is towards higher coordinates on the minus strand
                lo = max(0, start - b)
                hi = end + a
                seq = genome[chrom][lo:hi].complement.reverse.seq
            else:
                # '+' and unstranded '*' are treated the same way
                lo = max(0, start - a)
                hi = end + b
                seq = genome[chrom][lo:hi].seq

            # minimum sequence length
            if len(seq) >= minlen:
                name = fields[3] + '|' + '|'.join([chrom, str(start), str(end)])
                outfa.write('>' + name + '|' + strand + '\n')
                outfa.write(seq + '\n')
                written += 1

    return written
