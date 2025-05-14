#!/usr/bin/env python3

import sys
import re
import csv
import pandas as pd
from pandas.errors import EmptyDataError

def normalize_sample_name(name, prefix="CA"):
    # Dynamically build the regex based on the prefix
    pattern = rf'{re.escape(prefix)}-(.*?)(\.|$)'
    match = re.search(pattern, name)
    if not match:
        return None
    segment = match.group(1)
    # Remove hyphen followed by digit
    normalized = re.sub(r'-\d+', '', segment)
    return normalized

def update_matrix(input_file, matrix_file, string_species, output_file=None, prefix=""):
    # Read input values
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        raise ValueError("Input file must contain at least a header and one data line")

    headers = lines[0].strip().split()
    values = None
    
    # Try to find a line starting with 'average'
    for line in lines[1:]:
        parts = line.strip().split()
        if parts and parts[0].lower() == "average":
            values = list(map(float, parts[1:]))
            break

    # Fallback if no 'average' line is found
    if values is None:
        if len(lines) == 2:
            parts = lines[1].strip().split()
            try:
                values = list(map(float, parts[1:]))
            except ValueError:
                raise ValueError("Second line does not contain numeric values")
        else:
            raise ValueError("No line starting with 'average' found, and multiple data lines exist")

    # Normalize sample names
    sample_names = [normalize_sample_name(h, prefix) for h in headers]
    if None in sample_names:
        raise ValueError("Could not extract all sample names")

    # Load or create the matrix file
    try:
        df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    except (FileNotFoundError, EmptyDataError):
        df = pd.DataFrame()

    
    # Ensure row indices are strings
    df.index = df.index.astype(str)

    # Ensure species column exists
    if string_species not in df.columns:
        df[string_species] = pd.Series(dtype=float)

    # Add or update values
    for name, value in zip(sample_names, values):
        if name not in df.index:
            # Safely add a new row with NaNs for all columns
            df.loc[name] = pd.Series({col: float('nan') for col in df.columns})
        # Update the value
        df.at[name, string_species] = value

    # Drop any duplicated rows (just in case)
    df = df[~df.index.duplicated(keep='last')]

    # Output
    output_path = output_file if output_file else matrix_file
    df.to_csv(output_path, sep='\t')


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: update_matrix.py <input_file> <matrix_file> <string_species> [output_file]")
        sys.exit(1)

    input_file = sys.argv[1]
    matrix_file = sys.argv[2]
    species = sys.argv[3]
    output_file = sys.argv[4] if len(sys.argv) > 4 else None
    prefix = sys.argv[5] if len(sys.argv) > 5 else "CA"

    print(input_file)
    print(matrix_file)
    update_matrix(input_file, matrix_file, species, output_file, prefix)
