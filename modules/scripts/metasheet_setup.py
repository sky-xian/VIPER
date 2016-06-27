#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import pandas as pd

def updateMeta(config):
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',')
    config["comparisons"] = [c[5:] for c in metadata.columns if c.startswith("comp_")]
    config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'comp']
    config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
    config["ordered_sample_list"] = metadata.index
    return config
