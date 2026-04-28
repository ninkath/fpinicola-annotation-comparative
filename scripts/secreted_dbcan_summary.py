#!/usr/bin/env python3

import re
import pandas as pd
from collections import defaultdict

in_files = {
    "F_pinicola": snakemake.input.pinicola,
    "F_schrenkii": snakemake.input.schrenkii,
    "F_rosea": snakemake.input.rosea
}
outfile = snakemake.output.summary

fam_pattern = re.compile(r"(GH|GT|AA|CE|CBM|PL)[0-9]+(?:_[0-9]+)?")
preferred_cols = ["Recommend Results", "DIAMOND", "dbCAN_hmm", "dbCAN_sub"]

family_counts = defaultdict(dict)

for species, path in in_files.items():
    df = pd.read_csv(path, sep="\t")
    # No need for #ofTools filter here - already filtered + signal peptide from merge_secreted_cazymes.py
    
    present_cols = [c for c in preferred_cols if c in df.columns]
    if not present_cols:
        raise ValueError(f"No expected CAZy columns found in {path}")
    
    counts = defaultdict(int)
    
    for _, row in df.iterrows():
        found = set()
        for col in present_cols:
            val = row[col]
            if pd.isna(val) or str(val).strip() in ["", "-"]:
                continue
            for m in fam_pattern.finditer(str(val)):
                found.add(m.group(0))
        for fam in found:
            counts[fam] += 1
    
    for fam, n in counts.items():
        family_counts[fam][species] = n

out = pd.DataFrame.from_dict(family_counts, orient="index").fillna(0).astype(int)
out.index.name = "Family"
out = out.reset_index()

for sp in ["F_pinicola", "F_schrenkii", "F_rosea"]:
    if sp not in out.columns:
        out[sp] = 0

out = out[["Family", "F_pinicola", "F_schrenkii", "F_rosea"]]
out.to_csv(outfile, sep='\t', index=False)
