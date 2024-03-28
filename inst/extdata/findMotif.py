import re

def findMotif(peaksfa,motif='[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]'):
    # save in list
    relpos = []
    # open fa
    with open(peaksfa,'r') as seq:
        for line in seq:
            line = line.replace('\n','')
            if line.startswith('>'):
                next
            # search all motif
            pattern = re.compile(motif,flags=re.IGNORECASE)
            pos_res = pattern.findall(line)
            # save in list
            relpos.extend(pos_res)
    return relpos