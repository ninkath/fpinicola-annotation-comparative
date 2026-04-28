#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv(snakemake.input.presence, sep="\t")

species = ["F_pinicola", "F_schrenkii", "F_rosea"]

def label_row(row):
    present = [sp for sp in species if row[sp] == 1]
    return "+".join(present)

df["Combination"] = df.apply(label_row, axis=1)

summary = (
    df.groupby("Combination")
      .size()
      .reset_index(name="Count")
      .sort_values("Count", ascending=False)
)

summary.to_csv(snakemake.output.summary, sep="\t", index=False)
