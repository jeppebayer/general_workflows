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


# ────────────────────────────────────────────────────────────
# 1.  name normalisation
# ────────────────────────────────────────────────────────────
REMOVE_PATTERNS = re.compile(r'(?:^|[_-])(CA|K|GR)(?:[_-])+', flags=re.I)
TRAILING_HYPHEN_NUM = re.compile(r'-\d+')

def normalise(name: str, prefix: str = "CA") -> str:
    """Return a normalised sample name according to the rules."""
    # keep only the part after the first prefix match
    m = re.search(rf'{re.escape(prefix)}[_-](.*)', name, flags=re.I)
    seg = m.group(1) if m else name

    # drop duplicated prefixes like CA_, CA-  (case-insensitive)
    seg = REMOVE_PATTERNS.sub('', seg)

    # collapse ST-6J  →  STJ   (second hyphen followed by number)
    seg = TRAILING_HYPHEN_NUM.sub('', seg)

    # build final name  (no leading underscore for GR)
    if prefix.upper() == "GR":
        base = seg
    else:
        base = f"{prefix}_{seg}"

    return base


def deduplicate(names):
    """If duplicates occur, append last -C55 etc. to make them unique."""
    seen = {}
    out  = []
    for raw, base in names:
        if base not in seen:
            seen[base] = 1
            out.append(base)
        else:
            # restore the final -C55 (or similar) so it becomes unique
            tail = raw.split('-')[-1]
            unique = f"{base}-{tail}"
            while unique in seen:            # extreme corner-case
                unique += "_dup"
            seen[unique] = 1
            out.append(unique)
    return out


# ────────────────────────────────────────────────────────────
def update_matrix(input_path, matrix_path, species, output_path, prefix="CA"):
    # ---------- read input single-sample file ----------
    with open(input_path, encoding='utf-8') as f:
        header_line = f.readline().rstrip('\n')
        value_line  = f.readline().rstrip('\n')

    headers = header_line.split('\t')
    values  = list(map(float, value_line.split('\t')))
    if len(headers) != len(values):
        raise ValueError("Header and value counts differ")

    # ---------- normalise & deduplicate names ----------
    bases = [normalise(h, prefix) for h in headers]
    norm_names = deduplicate(list(zip(headers, bases)))

    # ---------- load matrix (or create empty) ----------
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
    df.to_csv(output_path, sep='\t')


# ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    in_file, matrix_file, species, out_file = sys.argv[1:5]
    pref = sys.argv[5] if len(sys.argv) > 5 else "CA"

    update_matrix(Path(in_file),
                  Path(matrix_file),
                  species,
                  Path(out_file),
                  pref)
