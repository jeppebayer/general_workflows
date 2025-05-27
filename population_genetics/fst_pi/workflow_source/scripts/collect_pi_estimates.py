#!/usr/bin/env python

import sys
import pandas as pd

#if len(sys.argv) < 5:
 #   print("Usage: python script.py <collection_file> <row_name> <column_name> <estimate>")
  #  #sys.exit(1)  # Exit with error

collection_file = sys.argv[1]
row_name = sys.argv[2]
column_name = sys.argv[3]
estimate = sys.argv[4]


# function for updating a cell with pi values
def update_cell(file_path, row_name, col_name, estimate, delimiter="\t"):
    """
    Updates a specific cell in a CSV/TSV file. If the row or column does not exist, it creates them.
    
    Parameters:
    - file_path: str, path to the CSV/TSV file.
    - row_name: str, row index name to search for.
    - col_name: str, column name to search for.
    - estimate: any, value to insert into the cell.
    - delimiter: str, delimiter used in the file (default: ',').
    """
    if not os.path.exists(file_path):
        print(f"File not found, creating new one.: {file_path}")
        df = pd.DataFrame()
    else:
        df = pd.read_csv(file_path, index_col=0, sep=delimiter)

    #try:
     #   # Load the file
      #  df = pd.read_csv(file_path, index_col=0, sep=delimiter)
       # print(df.columns.tolist())
    #except FileNotFoundError:
     #   # If the file does not exist, create a new DataFrame
      #  df = pd.DataFrame()

    # Ensure the DataFrame has at least one column
    if df.empty:
        df = pd.DataFrame(columns=[column_name])

    # Ensure the row exists
    if row_name not in df.index:
        df.loc[row_name] = [None] * len(df.columns)  # Add a new row

    # Ensure the column exists
    if col_name not in df.columns:
        df[col_name] = None  # Add a new column with NaN values

    # Update the cell with the estimate
    df.at[row_name, col_name] = float(estimate)

    # Save the modified DataFrame back to file
    df.to_csv(file_path, sep=delimiter)
	# thank you chat gbt


# run actual insert

update_cell(file_path=collection_file, row_name=row_name, col_name=column_name, estimate=estimate)

exit()