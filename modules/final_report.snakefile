#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#---------------------------------

from modules.scripts.viper_report import get_sphinx_report
from snakemake.utils import report
from modules.scripts.utils import getTargetInfo
from modules.scripts.utils import _copyMetaFiles

rule generate_report:
    input:
        getTargetInfo(config),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/" + config["token"] + ".html"
    message: "Generating VIPER report"
    benchmark:
        "benchmarks/" + config["token"] + "/generate_report.txt"
    run:
        sphinx_str = get_sphinx_report(config)
        report(sphinx_str, output[0], metadata="Molecular Biology Core Facilities, DFCI", **{'Copyrights:':"./viper/mbcf.jpg"})

rule copy_meta_files:
    input:
        config_file = config["config_file"], 
        meta_file = config["metasheet"]
    output:
        config_file = _copyMetaFiles(config)[0], 
        meta_file = _copyMetaFiles(config)[1]
    message: "Saving config and metasheet into analysis/{config[token]}"
    shell:
        "cp {input.config_file} {output.config_file} && "
        "cp {input.meta_file} {output.meta_file}"         
