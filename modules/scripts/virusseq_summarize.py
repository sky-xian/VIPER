#!/usr/bin/env python
"""Script to summarize the virusseq results--outputs to stdout:
SampleID, TranscriptID, Counts, FPKM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f "
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="virusseq.table file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f=open(options.file)
    #BURN
    tmp = f.readline()
    #RELYING on the table file already being sorted!
    lastID = ''
    #print header
    print("\t".join(["SampleID","TranscriptID","Counts","FPKM"]))
    for l in f:
        tmp = l.strip().split(',')
        if tmp[0] != lastID:
            #NEW SET
            ct = 1
            lastID = tmp[0]
            #PRINT the whole row
            print("\t".join(tmp))
            continue

        if ct < 5:
            #PRINT JUST the new info
            tmp[0] = ''
            print("\t".join(tmp))
            ct += 1

if __name__=='__main__':
    main()


