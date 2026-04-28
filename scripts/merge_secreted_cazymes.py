#!/usr/bin/env python3

import sys
import pandas as pd

dbcan_file = sys.argv[1]
signalp_file = sys.argv[2]
output_file = sys.argv[3]

def get_core_id(gene_id):
    """
    Extract Protein ID (number) from JGI-style names,
    or keep the full name for FUN_0000X-T1 IDs.
    """
    gene_id = str(gene_id)
    if gene_id.startswith('jgi|'):
        # # Split on | and take the 3rd element (the numeric ID)
        return gene_id.split('|')[2]
    elif gene_id.startswith('jgi-'):
        # Split on - and take the 3rd element (the numeric ID)
        return gene_id.split('-')[2]
    else:
        # F. pinicola annotation IDs (e.g. FUN_000001-T1)
        return gene_id

# 1. Read SignalP and collect the core ID of every secreted protein
secreted_core_ids = set()
with open(signalp_file, 'r') as f:
    for line in f:
        if "CS pos:" in line:
            raw_id = line.split()[0]
            core_id = get_core_id(raw_id)
            secreted_core_ids.add(core_id)

# 2. Read dbCAN and apply the consensus filter (>= 2 of 3 methods)
df_dbcan = pd.read_csv(dbcan_file, sep='\t')
df_dbcan = df_dbcan[df_dbcan["#ofTools"] >= 2].reset_index(drop=True)
gene_col = df_dbcan.columns[0]

# 3. Cross-check by extracting the core ID from dbCAN entries as well
mask = df_dbcan[gene_col].apply(lambda x: get_core_id(x) in secreted_core_ids)
df_secreted_cazymes = df_dbcan[mask]

# 4. Save the result
df_secreted_cazymes.to_csv(output_file, sep='\t', index=False)

print(f"Found {len(df_secreted_cazymes)} secreted CAZymes for this species.")
