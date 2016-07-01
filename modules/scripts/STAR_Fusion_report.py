#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

# @AUTHOR: Mahesh Vangala
# @Email: vangalamaheshh@gmail.com
# @Date: May,9,2016

import argparse
import sys
import pandas as pd
import os.path
import re

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fusionOutFiles', required=True, nargs='*', help="Provide fusion output filename (*final.abridged). " + 
        "Multiple filenames can be given as space separated values.")
    args = parser.parse_args()
    return args

def printFusionInfo( file_list ):
    main_df = pd.DataFrame()
    for cur_file in file_list:
        df = pd.read_table(cur_file, index_col=False)
        sample = re.sub("\..*", "", os.path.basename(cur_file))
        df["Sample"] = sample
        df["TotalReads"] = df["JunctionReads"] + df["SpanningFrags"]
        df["FusionName"] = df["#fusion_name"]
        main_df = main_df.append(df[["Sample","FusionName","JunctionReads","SpanningFrags","TotalReads","Splice_type", "LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint"]])
    print(main_df.to_csv(header=True, index=False))

if __name__ == "__main__":
    args = parseArgs()
    printFusionInfo(args.fusionOutFiles)
