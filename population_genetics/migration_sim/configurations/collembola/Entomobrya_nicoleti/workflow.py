#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('../../../workflow_source/'))
from workflow_source import *

gwf = migration_simulation()


# creates 2dSFS and estimates migration rates based on those using FastSimCoal


# conda create -n migration_fsc python R gwf pyyaml r-ggplot2 r-viridis

# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/migration_sim/configurations/collembola/Entomobrya_nicoleti
# conda activate migration_fsc
# gwf run
# gwf status
# gwf status -f summary
# gwf logs pi_remodelling_file
# gwf logs --stderr pi_remodelling_file
#⨯ pi_calculation_all_positions failed
#⨯ pi_add_context cancelled


