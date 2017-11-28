#!/usr/bin/env python
"""Script to generate a list of potential viruses in virusseq-genes.results

Method: grep for chrM transcripts that have an TPM > 1.0
NOTE: I probably should check against a db
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f virusseq.transcripts.gtf"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="virusseq gtf file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    #header
    hdr = f.readline().strip()
    print(hdr)
    for l in f:
        tmp = l.strip().split('\t')
        #ONLY look at chrM transcripts--where the viruses are placed
        if tmp[7] == 'chrM' and float(tmp[5]) >= 1.0:
            print("\t".join(tmp))

if __name__=='__main__':
    main()


