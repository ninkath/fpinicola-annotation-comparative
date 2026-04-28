#!/usr/bin/env python3

import re
import pandas as pd

species_files = snakemake.params.species
out_file = snakemake.output.summary

classes = ["GH", "GT", "AA", "CE", "CBM", "PL"]

# regex for family hits like GH36, AA3_2, CBM50, GT2, PL1, CE16, etc.
fam_pattern = re.compile(r"(GH|GT|AA|CE|CBM|PL)[0-9]+(?:_[0-9]+)?")

# columns worth scanning in overview.tsv
preferred_cols = ["Recommend Results", "DIAMOND", "dbCAN_hmm", "dbCAN_sub"]

rows = []

for species, path in species_files.items():
    df = pd.read_csv(path, sep="\t")
    df = df[df["#ofTools"] >= 2].reset_index(drop=True)   # consensus filter

    present_cols = [c for c in preferred_cols if c in df.columns]
    if not present_cols:
        raise ValueError(f"No expected CAZy columns found in {path}. Columns: {list(df.columns)}")

    counts = {c: 0 for c in classes}

    for _, row in df.iterrows():
        found = set()

        for col in present_cols:
            val = row[col]
            if pd.isna(val) or str(val).strip() == "-" or str(val).strip() == "":
                continue

            matches = fam_pattern.findall(str(val))
            # findall with groups returns only the class names if pattern has capture groups
            # so instead re-scan using finditer
            for m in fam_pattern.finditer(str(val)):
                fam = m.group(0)
                for cls in classes:
                    if fam.startswith(cls):
                        found.add(cls)
                        break

        for cls in found:
            counts[cls] += 1

    for cls in classes:
        rows.append({
            "Class": cls,
            "Species": species,
            "Count": counts[cls]
        })

out = pd.DataFrame(rows)
wide = out.pivot(index="Class", columns="Species", values="Count").reset_index()
wide = wide[["Class", "F_pinicola", "F_schrenkii", "F_rosea"]]
wide.to_csv(out_file, sep="\t", index=False)
