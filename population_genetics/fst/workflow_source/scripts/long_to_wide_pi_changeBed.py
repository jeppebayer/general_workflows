#!/usr/bin/env python3

#import packages
import pandas as pd
import numpy as np
import sys


# Ensure correct number of arguments are provided
if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <arg1> <arg2> <arg3>")
    sys.exit(1)

# Import positional arguments
sorted_pi_file = sys.argv[1]
output_wide_pi = sys.argv[2]
output_bed_pi = sys.argv[3]

# Print the arguments
print(f"Argument 1, sorted pi file input: {sorted_pi_file}")
print(f"Argument 2, pi output: {output_wide_pi}")
print(f"Argument 3, bed output: {output_bed_pi}")

# starting actual script

# Read data from a file (assuming 'data.txt' is the file name)
df = pd.read_csv(sorted_pi_file, sep='\\s+')
#df.info()

# Pivot the data to wide format
wide_df = df.pivot(index=['chrom', 'chromStart', 'chromEnd'], columns='pop_name', values='pi').reset_index()
#df = []	# reducing memory impact of dataframe
del df	# reducing memory impact of dataframe

# reindex to remove nan column
wide_df = wide_df.reindex(index=wide_df.index.difference([np.nan]), columns=wide_df.columns.difference([np.nan]))
#wide_df.info()

# fill nas with 0
wide_df_fill = wide_df.fillna('0')
#wide_df_fill.info()

#wide_df = []	# reducing memory impact of dataframe
del wide_df	# reducing memory impact of dataframe

# sort columns
# separate dataset to sort on pop cols
wide_df_fill_first = wide_df_fill[['chrom','chromStart','chromEnd']]
wide_df_fill_second = wide_df_fill[wide_df_fill.columns.difference(['chrom','chromStart','chromEnd'])].sort_index(axis = 'columns', key=lambda col: col.str.lower())
#wide_df_fill_second.info()

# combine split sets
wide_df_fill = pd.concat([ wide_df_fill_first, wide_df_fill_second ], ignore_index=False, axis=1)	#reusing df name

del wide_df_fill_second 	# reducing memory impact of dataframe
#wide_df_fill_second = []	# reducing memory impact of dataframe

# Write output to a file or print it
#output_file = output_wide_pi
# output_wide_pi = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.pi'
wide_df_fill.to_csv(output_wide_pi, sep='\t', index=False)


## Make bed file


# make comma separated list for "name" and "score" columns in bed format
wide_df_fill_concatPi = wide_df_fill.apply(lambda row: ", ".join(row.iloc[3:].astype(str)), axis = 'columns')
wide_df_fill_concatPops = ", ".join(wide_df_fill.columns[3:].astype(str))

del wide_df_fill	# reducing memory impact of dataframe
#wide_df_fill = []	# reducing memory impact of dataframe

# make new pd dataframe of score and name columns
new_df = pd.DataFrame({
    'name': 'pi' * len(wide_df_fill_concatPi),  # Repeat header for all rows
    'score': wide_df_fill_concatPi})

#del wide_df_fill_concatPops 	# reducing memory impact of dataframe
del wide_df_fill_concatPi	# reducing memory impact of dataframe
#wide_df_fill_concatPops = []	# reducing memory impact of dataframe
#wide_df_fill_concatPi = []	# reducing memory impact of dataframe

    
wide_df_bed = pd.concat([ wide_df_fill_first, new_df ], ignore_index=False, axis=1)

del wide_df_fill_first
del new_df
# write to bedfile
#output_file = output_bed_pi
# output_bed_pi = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.bed'
wide_df_bed.to_csv(output_bed_pi, sep='\t', index = False, header = False)
# output population order for the pi list
wide_df_fill_concatPops.to_csv(output_bed_pi.replace(".bed", "_populationOrder.txt"), sep='\t', index = False, header = False)

print('Finished remodelling data')