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

rule generate_report:
    input:
        getTargetInfo(config),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/report.html"
    message: "Generating VIPER report"
    run:
        sphinx_str = get_sphinx_report(config["comparisons"])
        report(sphinx_str, output[0], metadata="Molecular Biology Core Facilities, DFCI", **{'Copyrights:':"./viper/mbcf.jpg"})


