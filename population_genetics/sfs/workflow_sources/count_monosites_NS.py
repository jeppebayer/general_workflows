import sys
import pandas as pd

def calculate_probabilities(input_file, output_file):
    # Read the input tab-separated file into a pandas DataFrame
    data = pd.read_csv(input_file, sep="\t")
    
    # Group by chromosome and calculate the sums for Missense Probability and Synonymous Probability
    grouped_data = data.groupby("Chromosome").agg(
        N_count=("Missense Probability", "sum"),
        S_count=("Synonymous Probability", "sum")
    ).reset_index()
    
    # Save the results to a tab-separated file
    grouped_data.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")

# Define the input and output file paths
input_file = sys.argv[1] #"input_file.tsv"  # Replace with the path to your input file
output_file = sys.argv[2] #"output_file.tsv"  # Replace with the path to your output file

# Run the function
calculate_probabilities(input_file, output_file)

