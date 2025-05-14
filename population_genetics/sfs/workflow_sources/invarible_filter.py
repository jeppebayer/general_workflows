import sys
import pandas as pd

# Function to extract DP and TYPE values from the Nearest_Variant_Info field
def extract_info(info):
    # Split the information string by semicolon
    info_parts = info.split(';')
    dp_value = None
    type_value = None
    # Loop through the parts to find DP and TYPE values
    for part in info_parts:
        if part.startswith("DP="):
            dp_value = int(part.split('=')[1])
        if part.startswith("TYPE="):
            type_value = part.split('=')[1].replace("('", "").replace("',)", "")
    return dp_value, type_value

# Function to apply filters and decide whether to keep a site
def apply_filters(query_position, nearest_variant_position, dp_value, type_value, dp_lower, dp_upper):
    # Calculate the absolute distance
    absolute_distance = abs(query_position - nearest_variant_position)
    
    # Filter condition 1: absolute_distance should be <= 50
    if absolute_distance > 50:
        return False
    
    # Filter condition 2: DP should be within the given thresholds
    if dp_value < dp_lower or dp_value > dp_upper:
        return False
    
    # Filter condition 3: If type_value contains "del" or "ins" and absolute_distance <= 5, remove the site
    if ("del" in type_value or "ins" in type_value) and absolute_distance <= 5:
        return False
    
    return True

# Function to read depth thresholds from the provided file
def read_depth_thresholds(file_path):
    df = pd.read_csv(file_path, sep='\t')
    min_coverage_threshold = df['min_coverage_threshold'].values[0]
    max_coverage_threshold = df['max_coverage_threshold'].values[0]
    return min_coverage_threshold, max_coverage_threshold

# Function to process the TSV file line by line, apply filters, and generate a summary
def process_tsv_with_summary(variant_file, depth_threshold_file, output_file, summary_file):
    # Read DP thresholds from the depth threshold file
    dp_lower, dp_upper = read_depth_thresholds(depth_threshold_file)
    
    # Initialize a counter for the number of sites passing the filter
    passing_sites_count = 0

    # Open the input TSV file and the output file
    with open(variant_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write the header to the output file
        outfile.write("Query_Chrom\tQuery_Position\tDP\tTYPE\tAbsolute_Distance\n")

        # Skip the header line in the input file
        header = infile.readline()

        # Process each line one by one
        for line in infile:
            # Split the line into fields (tab-separated)
            fields = line.strip().split('\t')

            # Extract the relevant fields
            query_chrom = fields[0]
            query_position = int(fields[1])
            nearest_variant_position = int(fields[2])
            nearest_variant_info = fields[5]

            # Extract DP and TYPE from the Nearest_Variant_Info field
            dp_value, type_value = extract_info(nearest_variant_info)

            # Apply the filters
            if apply_filters(query_position, nearest_variant_position, dp_value, type_value, dp_lower, dp_upper):
                # If the site passes the filters, calculate the absolute distance and write the result
                absolute_distance = abs(query_position - nearest_variant_position)
                outfile.write(f"{query_chrom}\t{query_position}\t{dp_value}\t{type_value}\t{absolute_distance}\n")
                # Increment the counter for passing sites
                passing_sites_count += 1
    
    # Write the summary to the summary file
    with open(summary_file, 'w') as summary_outfile:
        summary_outfile.write(f"Total number of sites passing the filters: {passing_sites_count}\n")

# Example usage
variant_file = sys.argv[1]#'input_variants.tsv'  # Replace with your actual variants file path
depth_threshold_file = sys.argv[2]#'depth_thresholds.tsv'  # Replace with your actual depth threshold file path
output_file = sys.argv[3]#'filtered_output.tsv'  # Output file to store the filtered results
summary_file = sys.argv[4]#'summary_output.txt'  # Summary file to store the number of passing sites

# Process the file with filters and generate a summary
process_tsv_with_summary(variant_file, depth_threshold_file, output_file, summary_file)

