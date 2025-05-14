import sys
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <folded_output_file> <unfolded_output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    folded_output_file = sys.argv[2]
    unfolded_output_file = sys.argv[3]

    # Load the input file
    try:
        data = pd.read_csv(input_file, sep="\t")
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Summarize data by count_type and count_number
    try:
        folded_sfs = data[data['count_type'] == 'fold_count'].groupby('count_number')['number_of_sites'].sum().reset_index()
        unfolded_sfs = data[data['count_type'] == 'alt_count'].groupby('count_number')['number_of_sites'].sum().reset_index()
    except KeyError as e:
        print(f"Missing expected column in input file: {e}")
        sys.exit(1)

    # Write the results to output files
    try:
        folded_sfs.to_csv(folded_output_file, sep="\t", index=False, header=["count_number", "number_of_sites"])
        unfolded_sfs.to_csv(unfolded_output_file, sep="\t", index=False, header=["count_number", "number_of_sites"])
    except Exception as e:
        print(f"Error writing output files: {e}")
        sys.exit(1)

    print(f"Folded SFS written to: {folded_output_file}")
    print(f"Unfolded SFS written to: {unfolded_output_file}")

if __name__ == "__main__":
    main()

