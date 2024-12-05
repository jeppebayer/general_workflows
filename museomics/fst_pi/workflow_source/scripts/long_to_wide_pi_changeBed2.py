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
#!/usr/bin/env python3

# Define column types for more efficient memory usage
dtype_dict = {
    'chrom': 'category',
    'chromStart': 'int32',
    'chromEnd': 'int32',
    'pop_name': 'category',
    'pi': 'float32'  # Assuming 'pi' is a float value
}

# Read the input CSV in chunks to handle large files
chunk_size = 150  # Adjust chunk size as necessary
chunks = pd.read_csv(sorted_pi_file, sep=r'\s+', dtype=dtype_dict, chunksize=chunk_size)

# Initialize an empty list to hold pivoted dataframes
wide_dfs = []

for chunk in chunks:
    # Pivot the chunk
    pivot_chunk = chunk.pivot(index=['chrom', 'chromStart', 'chromEnd'], columns='pop_name', values='pi')
    wide_dfs.append(pivot_chunk)

# Concatenate all pivoted chunks
wide_df = pd.concat(wide_dfs).reset_index()

# Free up memory used by the chunks
del wide_dfs

# Fill NaN values with 0 and reorder columns
wide_df.fillna(0, inplace=True)

# Split data into two parts: the metadata columns and the population columns
metadata_columns = ['chrom', 'chromStart', 'chromEnd']
pop_columns = [col for col in wide_df.columns if col not in metadata_columns]
pop_columns_sorted = sorted(pop_columns, key=str.lower)

# Combine metadata and sorted population columns
wide_df = wide_df[metadata_columns + pop_columns_sorted]

# Write the wide dataframe to the output file in chunks to save memory
wide_df.to_csv(output_wide_pi, sep='\t', index=False)

# Generate BED file
# Create the "name" and "score" columns for BED format
wide_df['name'] = 'pi'
wide_df['score'] = wide_df[pop_columns_sorted].astype(str).apply(lambda row: ','.join(row), axis=1)

# Select only the relevant columns for BED format
bed_df = wide_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score']]

# Write BED output to file
bed_df.to_csv(output_bed_pi, sep='\t', index=False, header=False)

# Write population order to a separate file
with open(output_bed_pi.replace(".bed", "_populationOrder.txt"), 'w') as pop_order_file:
    pop_order_file.write(','.join(pop_columns_sorted))

print('Finished remodelling data')
