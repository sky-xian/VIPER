#!/usr/bin/env python
"""Script to generate a list of potential viruses in virusseq-transcript.gtf

Method: grep for chrM transcripts that have an FPKM > 0.0
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
    for l in f:
        tmp = l.strip().split('\t')
        #ONLY look at chrM transcripts--where the viruses are placed
        if tmp[0] == 'chrM':
            hit = False
            for a in tmp[-1].split(';'):
                foo = a.strip().split()
                #CHECK to see that FPKM > 0.0
                if foo and foo[0] == 'FPKM' and float(eval(foo[1])) > 0.0:
                    hit = True
                    break
            if hit:
                print("\t".join(tmp))

if __name__=='__main__':
    main()


