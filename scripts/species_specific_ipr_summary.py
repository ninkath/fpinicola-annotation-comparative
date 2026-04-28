#!/usr/bin/env python3

import pandas as pd

genes = pd.read_csv(snakemake.input.genes, sep="\t")

# -------------------------
# F. pinicola IPR
# -------------------------
pin_ipr = pd.read_csv(
    snakemake.input.pin_ipr,
    sep="\t",
    header=None,
    comment="#",
    usecols=[0, 11, 12],
    names=["Gene", "iprId", "iprDesc"]
)
pin_ipr = pin_ipr.dropna(subset=["iprId", "iprDesc"]).drop_duplicates()

# -------------------------
# JGI IPR tables
# -------------------------
def read_jgi_ipr(path, prefix):
    df = pd.read_csv(path, sep="\t", header=0)

    # first column is typically named "#proteinId"
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "proteinId"})

    df["proteinId"] = df["proteinId"].astype(str)
    df["Gene"] = df["proteinId"].apply(lambda x: f"jgi|{prefix}|{x}|")

    return (
        df[["Gene", "iprId", "iprDesc"]]
        .dropna(subset=["iprId", "iprDesc"])
        .drop_duplicates()
    )

sch_ipr = read_jgi_ipr(snakemake.input.sch_ipr, "Fompi3")
ros_ipr = read_jgi_ipr(snakemake.input.ros_ipr, "Fomro1")

all_rows = []

# exact match for pinicola
pin_genes = genes[genes["Species"] == "F_pinicola"].merge(
    pin_ipr, on="Gene", how="left"
)
all_rows.append(pin_genes)

# JGI species: match by proteinId parsed from OrthoFinder gene string
for sp, ipr_df, prefix in [
    ("F_schrenkii", sch_ipr, "Fompi3"),
    ("F_rosea", ros_ipr, "Fomro1"),
]:
    sub = genes[genes["Species"] == sp].copy()
    sub["proteinId"] = sub["Gene"].str.extract(
        rf"jgi\|{prefix}\|([0-9]+)\|"
    )[0].astype(str)

    ipr_df2 = ipr_df.copy()
    ipr_df2["proteinId"] = ipr_df2["Gene"].str.extract(
        rf"jgi\|{prefix}\|([0-9]+)\|"
    )[0].astype(str)

    merged = sub.merge(
        ipr_df2[["proteinId", "iprId", "iprDesc"]],
        on="proteinId",
        how="left"
    )
    all_rows.append(merged.drop(columns=["proteinId"]))

all_ipr = pd.concat(all_rows, ignore_index=True)

# clean descriptions
all_ipr["iprDesc"] = all_ipr["iprDesc"].astype(str).str.strip()
all_ipr = all_ipr[
    ~all_ipr["iprDesc"].isin(["-", r"\N", "", "nan", "None"])
]

# count unique genes per annotation
summary = (
    all_ipr.dropna(subset=["iprDesc"])
    .drop_duplicates(subset=["Species", "Gene", "iprDesc"])
    .groupby(["Species", "iprDesc"])["Gene"]
    .nunique()
    .reset_index(name="UniqueGenes")
    .sort_values(["Species", "UniqueGenes"], ascending=[True, False])
)

summary.to_csv(snakemake.output.table, sep="\t", index=False)
