#!/usr/bin/env python3

import sys
import re
import csv
import pandas as pd
from pandas.errors import EmptyDataError

def extract_segment(name, prefix="CA"):
    # Match prefix followed by either - or _
    pattern = rf'{re.escape(prefix)}[-_](.*?)(\.|$)'
    match = re.search(pattern, name)
    return match.group(1) if match else None

def normalize_sample_names(headers, prefix="CA"):
    segments = []
    short_names = []
    short_name_counts = {}

    # Step 1: Extract segments and short names
    for h in headers:
        segment = extract_segment(h, prefix)
        if segment is None:
            raise ValueError(f"Could not extract sample segment from: '{h}'")
        short_name = re.sub(r'-\w*\d+$', '', segment)  # Remove final -C123, -123 etc.
        segments.append(segment)
        short_names.append(short_name)
        short_name_counts[short_name] = short_name_counts.get(short_name, 0) + 1

    # Step 2: Apply logic
    normalized_names = []
    for segment, short_name in zip(segments, short_names):
        if short_name_counts[short_name] == 1:
            normalized_names.append(short_name)
        else:
            normalized_names.append(segment)
    return normalized_names


def update_matrix(input_file, matrix_file, string_species, output_file=None, prefix=""):
    # Read input values
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        raise ValueError("Input file must contain at least a header and one data line")

    #headers = lines[0].strip().split()
    headers = lines[0].rstrip('\n').split('\t')
    print("Headers:", headers)
    values = None
    
    # Try to find a line starting with 'average'
    for line in lines[1:]:
        #parts = line.strip().split()
        parts = line.rstrip('\n').split('\t')
        if parts and parts[0].lower() == "average":
            values = list(map(float, parts[1:]))
            break

    # Fallback if no 'average' line is found
    if values is None:
        if len(lines) == 2:
            #parts = lines[1].strip().split()
            parts = lines[1].rstrip('\n').split('\t')
            try:
                values = list(map(float, parts[1:]))
            except ValueError:
                raise ValueError("Second line does not contain numeric values")
        else:
            raise ValueError("No line starting with 'average' found, and multiple data lines exist")

    # Normalize sample names
    #sample_names = [normalize_sample_name(h, prefix) for h in headers]
    sample_names = normalize_sample_names(headers, prefix)
    print("Sample names:", sample_names)
    print("prefix:", prefix)
    if None in sample_names:
        for h in headers:
            norm = normalize_sample_names(h, prefix)
            if norm is None:
                print(f"⚠️ Failed to normalize: '{h}'")
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

    # Sort rows (index) and columns alphabetically
    df = df.sort_index().sort_index(axis=1)
    
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
