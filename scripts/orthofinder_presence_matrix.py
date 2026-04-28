#!/usr/bin/env python3

import pandas as pd

gene_counts = pd.read_csv(snakemake.input.gene_counts, sep="\t")

species = snakemake.params.species

# keep only orthogroup + species columns
cols = ["Orthogroup"] + species
df = gene_counts[cols].copy()

# convert counts to presence/absence
for sp in species:
    df[sp] = (df[sp] > 0).astype(int)

df.to_csv(snakemake.output.presence, sep="\t", index=False)
