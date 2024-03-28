import re

def CalculatePeaksRelativeDist(peaksfa,motif='[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]',mid=3):
    # save in list
    relpos = []
    # open fa
    with open(peaksfa,'r') as seq:
        for line in seq:
            line = line.replace('\n','')
            if line.startswith('>'):
                next
            peakmid = round(len(line)/2,ndigits=0)
            # search all motif
            pattern = re.compile(motif,flags=re.IGNORECASE)
            pos_res = pattern.finditer(line)
            # shift to A site position
            allpos = [pos.start() + mid for pos in pos_res]
            # calculate relative to peak center bases
            relposition = [pos - peakmid for pos in allpos]
            # save in list
            relpos.extend(relposition)
    return relpos