#!/bin/env python3
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('/faststorage/project/EcoGenetics/general_workflows/population_genetics/sfs/workflow_sources'))
from workflow_sources_multirun import *
configs = glob.glob('./*config.y*ml')
gwf = Workflow()
gwf = sfs_sp_build(configs[0],gwf)
for config in configs:
    gwf = sfs_pop_workflow(config_file = config,gwf = gwf)
    
