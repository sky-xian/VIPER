#!/usr/bin/env python
"""Script to collect the FPKM results (isoforms.fpkm_tracking, and maybe 
later genes.fpkm_tracking) across all samples
Outputs to stdout:
(Transcript/Gene)ID, Sample1, ..., SampleN
XXX_TRANSCRIPT_ID, SAMPLE1_FPKM, ..., SampleN_FPKM
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE1] -f [FPKM FILE2] ... -f [FPKM FILE N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of isoform files")
    optparser.add_option("-i", "--id_col", default=0, help="Column in the files that correspond to the id (default: 0)")
    optparser.add_option("-c", "--tpm_col", default=5, help="Column in the files that correspond to the TPM (default: 5)")
    optparser.add_option("-n", "--name_col", default=1, help="Column in the files that correspond to the gene_name")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    #rename:
    iid_col = int(options.id_col)
    tpm_col = int(options.tpm_col)
    name_col = int(options.name_col)

    #build up matrix 
    matrix = {}
    file_f = options.files[0]
    #INIT the samples list w/ the first sample name
    samples = [file_f.split("/")[-1].split(".")[0]]
    file_f = open(file_f)
    #BURN: read the header
    l = file_f.readline()
    for l in file_f:
        tmp = l.strip().split("\t")
        gene_name = tmp[name_col]
        #NOTE: these are duplicated--e.g. NOTCH1_NOTCH1--remove duplicated name
        gene_name = gene_name.split("_")[0]
        #initialize: matrix to {ID: (gene_name, [tpm_1, 0, 0, ..., 0])}
        matrix[tmp[iid_col]] = (gene_name, 
                                [tmp[tpm_col]] + [0]*(len(options.files)-1))
    file_f.close()

    #DEAL with rest:
    for (i, file_f) in enumerate(options.files[1:]):
        #APPEND to the samples name
        samples.append(file_f.split("/")[-1].split(".")[0])

        file_f = open(file_f)
        #BURN: read the header
        l = file_f.readline()
        for l in file_f:
            tmp = l.strip().split("\t")
            if (tmp[iid_col] in matrix):
                foo = matrix[tmp[iid_col]]
                #add to the tpm list 
                foo[1][i+1] = tmp[tpm_col]
                #ADD value
                matrix[tmp[iid_col]] = foo
            else:
                #NEW element--does this ever happen??
                foo = [0]*len(options.files)
                foo[i+1] = tmp[tpm_col]
                gene_name = tmp[name_col]
                #NOTE: remove duplicate name
                gene_name = gene_name.split("_")[0]
                matrix[tmp[iid]] = (gene_name, foo)
        file_f.close()

    #OUTPUT:
    print(",".join(["Transcript_ID"] + samples))

    for k in matrix.keys():
        gene_name = matrix[k][0]
        tpms = matrix[k][1]
        print(",".join([gene_name]+tpms))

if __name__=='__main__':
    main()


