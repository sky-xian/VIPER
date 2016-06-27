#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import yaml
import os, sys, subprocess

def updateConfig(config):
    ref_info = _getRefInfo()
    for k,v in ref_info.items():
        config[k] = v
    config["samples"] = config["the_samples"]
    config["config_file"] = "config.yaml" # trick to force rules on config change
	
    for k in ["RPKM_threshold","min_num_samples_expressing_at_threshold", 
                "SSnumgenes","SFnumgenes","num_kmeans_clust","filter_mirna","snp_scan_genome"]:
        config[k] = str(config[k])

    config = _addExecPaths(config)
    return config

def _getRefInfo():
    with open("ref.yaml","r") as ref_file:
        ref_info = yaml.safe_load(ref_file)
    return ref_info


def _addExecPaths(config):
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    #NEED the following when invoking python2 (to set proper PYTHONPATH)
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'viper_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = conda_path + '/python-2.7.9-3/bin/python2.7'

    if not "rseqc_path" in config or not config["rseqc_path"]:
        config["rseqc_path"] = conda_path + '/rseqc-2.6.2-0/bin'

    if not "picard_path" in config or not config["picard_path"]:
        config["picard_path"] = 'picard'

    if not "varscan_path" in config or not config["varscan_path"]:
        config["varscan_path"] = 'varscan'
    
    return config

