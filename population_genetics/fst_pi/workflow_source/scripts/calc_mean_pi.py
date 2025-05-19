#!/usr/bin/env python

import sys
import pandas as pd
from workflow_templates import *

#if len(sys.argv) < 5:
 #   print("Usage: python script.py <collection_file> <row_name> <column_name> <estimate>")
  #  #sys.exit(1)  # Exit with error

greenedalf_file = sys.argv[1]
species_short = sys.argv[2]
landcover_type = sys.argv[3]
ouput_dir = sys.argv[4]



pi_mean_greened(greenedalf_file, species_short, landcover_type, ouput_dir)

