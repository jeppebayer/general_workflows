#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../workflow_source/'))
from workflow_source import *

###########
###
####    This wf will run tests on all fastq files in listed directories within the yaml file
##
########

gwf = fastq_test_wf()

# conda create -n fastq_test_env python gwf pyyaml fastq_utils
# cd /home/anneaa/EcoGenetics/general_workflows/fastq_test_gwf_pipeline/configurations
# conda activate fastq_test_env
## conda env export > fastq_test_environment.yaml
# gwf run
# gwf status
# gwf status -f summary
# gwf logs --stderr make_genome_fai



