#!/usr/bin/env python

import sys
import pandas as pd
from workflow_templates import *

#if len(sys.argv) < 5:
 #   print("Usage: python script.py <collection_file> <row_name> <column_name> <estimate>")
  #  #sys.exit(1)  # Exit with error

def pi_mean_greened(input_file_greened: str, species_short: str, landscape_type: str, outplace:str) -> str:
    """Returns sum and mean of pi, wattersons theta and tajimas D, and exports them to each file with said name."""
    
    import pandas as pd
    # Load the input TSV
    df = pd.read_csv(input_file_greened, sep="\t")
    
    # Clean column headers: remove specific string from column names
    clean_str = "filtered_"  # ← replace this with the string you want to remove
    df.columns = [col.replace(clean_str, "") for col in df.columns]
    
    # Define the patterns to search for
    patterns = ["theta_pi", "theta_watterson", "tajimas_d"]
    
    # Process each pattern
    for pattern in patterns:
        # Select matching columns
        matching_cols = [col for col in df.columns if pattern in col]
        
        if not matching_cols:
            print(f"No columns matched for pattern: {pattern}")
            continue
        
        # Subset the DataFrame
        subset = df[matching_cols]
        
        # Compute sum and mean
        summed = subset.sum()
        averaged = subset.mean()
        counted = subset.count()
        
        # Combine into a new DataFrame
        result_df = pd.DataFrame([counted, summed, averaged], index=["count", "sum", "average"])
        
        # Write to output
        result_df.to_csv(f"{outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_summary.tsv", sep="\t")
        print(f"Summary written to {outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_summary.tsv, later changed to : {outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_mean.tsv")


greenedalf_file = sys.argv[1]
species_short = sys.argv[2]
landcover_type = sys.argv[3]
ouput_dir = sys.argv[4]



pi_mean_greened(greenedalf_file, species_short, landcover_type, ouput_dir)

exit()