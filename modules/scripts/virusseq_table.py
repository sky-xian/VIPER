#!/usr/bin/env python
"""Script to collect the virusseq results across all samples. outputs to stdout
SampleID, TranscriptID, Counts, TPM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [RSEM FILE_1] -f [RSEM FILE_2] ...-f [RSEM FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of virusseq.filtered.genes files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
    sampleIDs=[n.strip().split("/")[-1].split('.')[0] for n in options.files]

    print(",".join(["SampleID","TranscriptID","ExpectedCounts","TPM"]))

    for (i, rsem_f) in enumerate(options.files):
        #READ in tpm
        f = open(rsem_f)
        _dict = {}
        #burn the header
        l = f.readline()
        for l in f:
            tmp = l.strip().split("\t")

            transcript_id = tmp[0]
            expected_ct = tmp[4]
            tpm = tmp[5]
            _dict[transcript_id] = [expected_ct, tpm]
        #print(_dict)
        f.close()

        #COMPOSE TABLE
        table = sorted(_dict.items(), key=lambda x: x[1][0], reverse=True)

        #OUTPUT:
        for r in table:
            print(",".join([sampleIDs[i], r[0], str(r[1][0]), str(r[1][1])]))

if __name__=='__main__':
    main()


