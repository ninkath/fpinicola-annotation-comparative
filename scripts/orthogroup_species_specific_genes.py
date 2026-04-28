#!/usr/bin/env python3

import pandas as pd

orthogroups = pd.read_csv(snakemake.input.orthogroups, sep="\t")
presence = pd.read_csv(snakemake.input.presence, sep="\t")

species_cols = ["F_pinicola", "F_schrenkii", "F_rosea"]

merged = orthogroups.merge(
    presence,
    on="Orthogroup",
    suffixes=("_genes", "_present")
)

rows = []

for _, row in merged.iterrows():
    present = {sp: int(row[f"{sp}_present"]) for sp in species_cols}
    n_present = sum(present.values())

    if n_present == 1:
        species = [sp for sp in species_cols if present[sp] == 1][0]
        genes = row[f"{species}_genes"]

        if pd.isna(genes) or str(genes).strip() == "":
            continue

        for gene in [g.strip() for g in str(genes).split(",")]:
            if gene:
                rows.append({
                    "Orthogroup": row["Orthogroup"],
                    "Species": species,
                    "Gene": gene
                })

out = pd.DataFrame(rows)

if out.empty:
    out = pd.DataFrame(columns=["Orthogroup", "Species", "Gene"])

out.to_csv(snakemake.output.table, sep="\t", index=False)
