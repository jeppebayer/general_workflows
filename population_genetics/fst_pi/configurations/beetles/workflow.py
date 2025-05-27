#!/bin/env python3
import sys, os, glob
sys.path.insert(0, os.path.realpath('../../workflow_source/'))
from workflow_source import *
#os.environ['OPENBLAS_NUM_THREADS'] = '18'

# make list of ymls
configs = glob.glob('./*/*/*.config.y*ml')
configs = [file for file in configs if "deprec" not in file]

# if certain species should not be run:
#configs = [file for file in configs if "OrcCin" not in file]
#configs = [file for file in configs if "conservation_agriculture" not in file]
#configs = [file for file in configs if "EntNic" in file]
#configs = [file for file in configs if "PogFla" in file]
#configs = [file for file in configs if "OrcVil" in file]
#configs = [file for file in configs if "OrcCin" in file]
#configs = [file for file in configs if "conservation_agriculture" in file]
#configs = [file for file in configs if "conventional_agriculture" in file]
#configs = [file for file in configs if "grassland" in file]
print(configs)


# loop over them
gwf = Workflow()
for config in configs:
    gwf = fst_and_pi_wf(config_file = config, gwf = gwf)


# conda create -n ecogen_neutral_diversity_wf python samtools bamtools vcftools R gwf bedtools bcftools pyyaml r-ggplot2 r-viridis
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/fst_pi/configurations/collembola/Entomobrya_nicoleti/grassland
# conda activate ecogen_neutral_diversity_wf
# gwf run
# gwf status
# gwf status -f summary
# gwf status -s failed
# gwf logs pi_remodelling_file
# gwf logs --stderr pi_remodelling_file
#⨯ pi_calculation_all_positions failed
#⨯ pi_add_context cancelled


