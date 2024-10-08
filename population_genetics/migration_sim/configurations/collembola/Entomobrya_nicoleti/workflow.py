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
# gwf logs --stderr pair_2DSFS_map_target_519
#⨯ pi_calculation_all_positions failed
#⨯ pi_add_context cancelled

#36 min 2 cores 1 gb
#1t min 2 cores 1 gb
# similar with 1 core. huge time difference between jobs. some with 2 cores are still running 8h, while most finish quickly with both 1 or 2 coeres.
# could use -y parameter, indication that after X rounds of not inproving likelyhood, go back to the previous parameters improving likelyhood and retry.

