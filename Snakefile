from os import path

import numpy as np
import pandas as pd


################################################################################
# Globals                                                                      #
################################################################################

samples = (pd.read_csv(config['samples'], sep='\t')
             .set_index('sample', drop=False))

is_paired = config["options"]["paired"]
is_no_index = config["options"]["no_index"]

# default path for analysis
data_path = config["path"]["data"]
pipe_path = config["path"]["pipeline"]

# directory for output files
fastq_dir = path.join(data_path, 'raw')
fusion_dir = path.join(data_path, 'MapSplice')

qc_dir = path.join(data_path, 'qc')
log_dir = path.join(data_path, 'log')


################################################################################
# Functions                                                                    #
################################################################################

def get_samples():
    return set(samples['sample'])


def format_options(options):
    return ' '.join(options or [])

################################################################################
# Rules                                                                        #
################################################################################
rule all:
     input:
         expand(path.join(fusion_dir, '{sample}', 'fusions_well_annotated.txt'), 
         	sample = get_samples()),
         expand(path.join(fusion_dir, '{sample}', 'alignments.sort.bam'),
         	sample = get_samples())
         #path.join(qc_dir, 'multiqc_report.html')


include: path.join(pipe_path, 'rules/qc.smk')
include: path.join(pipe_path, 'rules/MapSplice.smk')






