from pyfaidx import Fasta

# define function
def GetPeaksFa(peak,genomefie,outfasta,type='bed',a=0,b=0,minlen=5):
    # load geneome
    genome = Fasta(genomefie)
    outfa = open(outfasta,'w')
    # process peaks file
    with open(peak,'r') as bed:
        for line in bed:
            fields = line.replace('\n','').split()
            chr = fields[0]
            strand = fields[5]
            start = int(fields[1])
            end = int(fields[2])
            if type == 'bed':
                end = end
            elif type == 'narrowpeak':
                end = end + 1
            else:
                print('please mark your peaks file type(bed/narrowpeak)')
                break
            # extend upstream and downstram peaks
            if strand == '+':
                seq = genome[chr][(start - a):(end + b)].seq
            elif strand == '-':
                seq = genome[chr][(start + a):(end -b)].complement.reverse.seq
            else:
                seq = genome[chr][(start - a):(end + b)].seq
            # mimimum seq length
            if len(seq) >= minlen:
                # fa name
                name = fields[3] + '|' + '|'.join([chr,str(start),str(end)])
                # output seq
                outfa.write('>' + name + '|' + strand + '\n')
                outfa.write(seq + '\n')

    print("Sequence has been written!")


    outfa.close()