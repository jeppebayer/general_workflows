#!/usr/bin/env python3
"""
update_matrix.py  <input_file>  <matrix_file>  <species>  <output_file>  [prefix]

* input_file   – single-line header + single-line values, TAB separated
* matrix_file  – existing TSV matrix (or new file)
* species      – name of the column to write/update (e.g. EntNic_CA)
* output_file  – path for the updated matrix
* prefix       – prefix used in normalisation (default "CA")
"""

import sys, re, pandas as pd
from pathlib import Path
from pandas.errors import EmptyDataError
#from collections import Counter
from collections import defaultdict
from filelock import FileLock
from pathlib import Path
import ast

# ────────────────────────────────────────────────────────────
# 1.  name normalisation
# ────────────────────────────────────────────────────────────
REMOVE_PATTERNS = re.compile(r'(?:^|[_-])(CA|K|GR)(?:[_-])+', flags=re.I)
TRAILING_HYPHEN_NUM = re.compile(r'-\d+')

def normalise(name: str, prefix: str = "CA") -> str:
    """Return a normalised sample name according to the rules."""
    # Extract everything after the first match of the prefix
    m = re.search(rf'{re.escape(prefix)}[_-](.*)', name, flags=re.I)
    seg = m.group(1) if m else name

    # Remove duplicated prefix-like segments: CA_, K-, GR_, etc.
    seg = REMOVE_PATTERNS.sub('', seg)

    # Remove all hyphen/underscore + trailing alphanum segments (e.g., -6J, _C55, -84)
    seg = re.sub(r'[-_][A-Za-z]*\d+', '', seg)
    seg = re.sub("ae", 'A', seg)
    seg = re.sub("aa", 'A', seg)
    seg = re.sub("oe", 'O', seg)

    # Final name
    if prefix.upper() == "GR":
        base = seg
    else:
        base = f"{prefix}_{seg}"

    return base




def deduplicate(names):
    """Ensure all normalized names are unique. If duplicates occur, all get suffixes."""
    base_to_raws = defaultdict(list)
    for raw, base in names:
        base_to_raws[base].append(raw)

    output = []
    for base, raws in base_to_raws.items():
        if len(raws) == 1:
            output.append(base)
        else:
            # Multiple entries → all get unique suffix
            for raw in raws:
                tail = raw.split('-')[-1]
                unique = f"{base}-{tail}"
                output.append(unique)

    return output

# ────────────────────────────────────────────────────────────
def update_matrix(input_path, matrix_path, species, output_path, prefix="CA", ignore_list = None):
    # ---------- read input single-sample file ----------
    print("inputfile is: ", input_path)

    with open(input_path, encoding='utf-8') as f:
        lines = f.readlines()
        
    header_line = lines[0].rstrip('\n')
    headers = header_line.split('\t')
        #header_line = f.readline().rstrip('\n')
        #value_line  = f.readline().rstrip('\n')

    headers = header_line.split('\t')
    first_data_line = lines[1].rstrip('\n').split('\t')
    #first_data_line = value_line.split('\t')
    #values_raw = value_line.split('\t')

    try:
        # Try converting all values
        #values = list(map(float, values_raw))
        values = list(map(float, first_data_line))
    except ValueError:
        # If conversion fails, assume first column is a row name and skip it
        headers = headers[1:]
        # Search for the correct row (row with average)
        target_rowname="average"
        #for line in f:
        for line in lines[1:]:
            parts = line.rstrip('\n').split('\t')
            if parts[0] == target_rowname:
                values = list(map(float, parts[1:]))
                break
        else:
            raise ValueError(f'Row name "{target_rowname}" not found in file.')
        
        #values = list(map(float, values_raw[1:]))

    #headers = header_line.split('\t')
    #values  = list(map(float, value_line.split('\t')))

    if len(headers) != len(values):
        raise ValueError("Header and value counts differ")

    # ---------- normalise & deduplicate names ----------
    bases = [normalise(h, prefix) for h in headers]
    norm_names = deduplicate(list(zip(headers, bases)))

    # -------- substutute species name followed by _ ---------
    sp_short = species.split('_')[0]
    norm_names = [re.sub(re.escape(f'{sp_short}_'), '', name) for name in norm_names]
    # -------- substutute data type name ---------
    norm_names = [re.sub(re.escape(".neutral.tajimas_d"), '', name) for name in norm_names]
    norm_names = [re.sub(re.escape(".neutral.theta_pi"), '', name) for name in norm_names]
    norm_names = [re.sub(re.escape(".neutral.theta_watterson"), '', name) for name in norm_names]

    # ------------ exclude pops from exclude list -----------
    if isinstance(ignore_list, list) and ignore_list:
        # Filter out excluded names
            # Remove None values and ensure all are strings
        ignore_list = [str(name) for name in ignore_list if name is not None]
        # now chage excluded names to match other names
        ignore_list = [re.sub(re.escape("oe"), "O", name) for name in ignore_list]
        ignore_list = [re.sub(re.escape("aa"), "A", name) for name in ignore_list]
        ignore_list = [re.sub(re.escape("ae"), "A", name) for name in ignore_list]
        filtered = [(n, v) for n, v in zip(norm_names, values) if n not in ignore_list]
        if filtered:
            norm_names, values = zip(*filtered)
        else:
            norm_names, values = [], []

    # ---------- load matrix (or create empty) ----------

    # ─────────────── 🔐 FILE LOCK SECTION ───────────────
    lock_path = str(matrix_path) + '.lock'
    lock = FileLock(lock_path)

    #lock_path = str(matrix_path) + '.lock'
    try:
        with lock: # FileLock(lock_path):

            try:
                df = pd.read_csv(matrix_path, sep='\t', index_col=0)
            except (FileNotFoundError, EmptyDataError):
                df = pd.DataFrame()

            # ensure column exists
            if species not in df.columns:
                df[species] = pd.NA

            # ensure all rows exist
            for rn in norm_names:
                if rn not in df.index:
                    df.loc[rn] = pd.NA

            # ---------- write values ----------
            for rn, val in zip(norm_names, values):
                df.at[rn, species] = val

            # ---------- sort rows + columns (case-insensitive) ----------
            df = df.reindex(sorted(df.index,   key=lambda s: s.lower()))        # rows
            df = df.reindex(sorted(df.columns, key=lambda s: s.lower()), axis=1)  # columns

            #df = df.reindex(sorted(df.index, key=lambda s: s.lower()))
            #cols = [c for c in df.columns if c != species]
            #cols_sorted = sorted(cols, key=lambda s: s.lower()) + [species]
            #df = df[cols_sorted]

            # ---------- save ----------
            #df.to_csv(output_path, sep='\t')
            df.to_csv(matrix_path, sep='\t')
    finally:
        #if FileLock(lock_path).is_locked:
            #FileLock(lock_path).release()
        if lock.is_locked:
            lock.release()
    
    #Path(out_file).touch()


# ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    
    in_file, matrix_file, species, out_file = sys.argv[1:5]
    pref = "CA"
    ignore_list = []

    # Optional 5th argument: prefix
    if len(sys.argv) >= 6:
        pref = sys.argv[5]

    # Optional 6th argument: ignore list as Python-style list string
    if len(sys.argv) >= 7:
        try:
            ignore_list = ast.literal_eval(sys.argv[6])
            if not isinstance(ignore_list, list):
                raise ValueError("Ignore list is not a list.")
        except Exception as e:
            print(f"Error parsing ignore list: {e}")
            sys.exit(1)

    update_matrix(Path(in_file),
              Path(matrix_file),
              species,
              Path(out_file),
              pref,
              ignore_list)
    
    exit()
