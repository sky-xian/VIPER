#!/usr/bin/env python
"""Script to collect the trust results across all samples. outputs to stdout
SampleID, AssemblyCount, TotalCount, CPK
Where:
AssemblyCount = # of lines in trust .fa (account for headers if any) / 2
**Don't forget to divide by 2!
Total count = 1st number in est_lib_size (4th field of first line)

CPK = (AssemblyCount*1000)/TotalCount
We graph CPK as a box plot and put into viper report

OUTPUTS to stdout
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [trust .fa output FILE_1] -f [trust .fa output FILE_2] ...-f [trust .fa output FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of trust .fa files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    #TRY to infer the SAMPLE NAMES--SAMPLE.sorted.bam.fa
    sampleIDs=[f.strip().split("/")[-1].split('.')[0] for f in options.files]
    #print(sampleIDs)
    
    print(",".join(["SampleID","AssemblyCount","TotalCount","CPK"]))

    for (i, trust_f) in enumerate(options.files):
        #get totalCount -- NOTE: assuming no header!!
        f = open(trust_f)
        tmp = f.readline().strip().split("+")
        #first num of 4th field
        totalCount = int(tmp[3].split("=")[1].split("-")[0])
        f.close()

        #calc: numLines
        #ref: https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
        num_lines = sum(1 for line in open(trust_f))
        assemblyCount = num_lines / 2.0
        
        CPK = (assemblyCount*1000)/totalCount
        print(",".join([sampleIDs[i], "%.1f" %assemblyCount, str(totalCount), "%.2f" % CPK]))

if __name__=='__main__':
    main()


