#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../workflow_source/'))
from workflow_source import *

gwf = fst_and_pi_wf()

# conda create -n ecogen_neutral_diversity_wf python samtools bamtools vcftools R gwf bedtools bcftools pyyaml
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/fst/configurations/collembola/Pogonognathellus_flavescens/
# conda activate ecogen_neutral_diversity_wf
# gwf run
# gwf status
# gwf status -f summary
# gwf logs --stderr make_genome_fai


