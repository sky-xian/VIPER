#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#------------------------------------

import yaml
import os, sys, subprocess

def updateConfig(config):
    ref_info = _getRefInfo()
    for k,v in ref_info.items():
        config[k] = v
    config["config_file"] = "config.yaml" # trick to force rules on config change
	
    for k in ["RPKM_threshold","min_num_samples_expressing_at_threshold", 
                "numgenes_plots","num_kmeans_clust","filter_mirna"]:
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
        #config["python2"] = conda_path + '/python-2.7.13-0/bin/python2.7'
        config["python2"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin', 'python2.7')

    if not "rseqc_path" in config or not config["rseqc_path"]:
        #config["rseqc_path"] = conda_path + '/rseqc-2.6.2-0/bin'
        config["rseqc_path"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin')

    if not "picard_path" in config or not config["picard_path"]:
        config["picard_path"] = 'picard'

    if not "varscan_path" in config or not config["varscan_path"]:
        config["varscan_path"] = 'varscan'
   
    # update RSEM executable only if rsem_ref key has value from user
    if "rsem_ref" in config and config["rsem_ref"]:
        config["rsem_path"] = conda_root + '/envs/rsem/bin'
        config["seurat_path"] = conda_root + '/envs/seurat/bin' 

    if "analysis_token" in config and config["analysis_token"]:
        config["token"] = config["analysis_token"]
    else:
        config["token"] = "summary_reports"
    
    return config

