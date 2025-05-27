#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../../workflow_source/'))
from workflow_source import *

config = glob.glob('*config.y*ml')[0]
gwf = Workflow()
gwf = fst_and_pi_wf(config_file = config, gwf=gwf)

# conda create -n ecogen_neutral_diversity_wf python samtools bamtools vcftools R gwf bedtools bcftools pyyaml r-ggplot2 r-viridis
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/fst_pi/configurations/collembola/Entomobrya_nicoleti/grassland
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/fst_pi/configurations/collembola
# conda activate ecogen_neutral_diversity_wf
# gwf run
# gwf status
# gwf status -f summary
# gwf status -s failed
# gwf logs pi_remodelling_file
# gwf logs --stderr pi_remodelling_file
#⨯ pi_calculation_all_positions failed
#⨯ pi_add_context cancelled


