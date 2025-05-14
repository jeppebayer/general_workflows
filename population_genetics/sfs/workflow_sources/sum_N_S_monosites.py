import sys
import csv

def sum_columns_by_name(file_path, column_names):
    # Initialize a dictionary to store sums for each column
    column_sums = {name: 0 for name in column_names}

    with open(file_path, 'r') as file:
        # Use csv.DictReader with '\t' as the delimiter for TSV files
        reader = csv.DictReader(file, delimiter='\t')
        
        # Iterate through each row in the file
        for row in reader:
            for name in column_names:
                # Safely convert to float and add to the sum
                column_sums[name] += float(row[name])
    
    return column_sums

# Specify the file path and columns to sum
file_path = sys.argv[1] #'data.tsv'
columns_to_sum = ['Missense Probability', 'Synonymous Probability']

# Calculate the sums
sums = sum_columns_by_name(file_path, columns_to_sum)

# Print the results
for column, total in sums.items():
    print(f"Sum of {column}: {total}")

