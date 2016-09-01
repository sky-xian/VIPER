#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------------

import pandas as pd

def updateMeta(config):
    _sanity_checks(config)
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')
    config["comparisons"] = [c[5:] for c in metadata.columns if c.startswith("comp_")]
    config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'comp']
    config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
    config["ordered_sample_list"] = metadata.index
    return config


def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', ' ':'_', '/':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()

