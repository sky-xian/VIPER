#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 09, 2016
#---------------------------

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule rsem_align:
    input:
        getFastq
    output:
        rsem_transprotected("analysis/RSEM/{sample}/{sample}.isoform.results")
        
