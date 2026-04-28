#!/usr/bin/env python3

import re
import pandas as pd
from collections import defaultdict

species_files = snakemake.params.species
out_file = snakemake.output.summary

# Match families like GH16, AA3_2, CBM50, CE4, GT2, PL1 etc.
fam_pattern = re.compile(r"(GH|GT|AA|CE|CBM|PL)[0-9]+(?:_[0-9]+)?")
preferred_cols = ["Recommend Results", "DIAMOND", "dbCAN_hmm", "dbCAN_sub"]

family_counts = defaultdict(dict)

for species, path in species_files.items():
    df = pd.read_csv(path, sep="\t")
    df = df[df["#ofTools"] >= 2].reset_index(drop=True)   # consensus filter
    present_cols = [c for c in preferred_cols if c in df.columns]

    if not present_cols:
        raise ValueError(f"No expected CAZy columns found in {path}. Columns: {list(df.columns)}")

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
out.to_csv(out_file, sep="\t", index=False)
