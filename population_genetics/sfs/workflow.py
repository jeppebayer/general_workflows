#!/bin/env python3
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('/faststorage/project/EcoGenetics/general_workflows/population_genetics/sfs/workflow_sources'))
from workflow_sources import *
configs = glob.glob('./configuratons/*config.y*ml')
for config in configs:
    gwf = purge_dup_for_assembly_workflow(config_file = config)
    
