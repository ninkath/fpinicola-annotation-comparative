#!/usr/bin/env python3
"""
Prepare plotting tables for a NUCmer pairwise dotplot using cumulative
genome coordinates.

Steps:
  1. Read NUCmer show-coords output (rclTH format).
  2. Filter alignments by minimum length and minimum identity.
  3. Reorder query sequences greedily so each reference sequence is paired
     with its best-matching query along the diagonal.
  4. Convert per-sequence coordinates to cumulative genome coordinates.
  5. Write two tables:
       - plot_coords.tsv: one row per retained alignment, with absolute
         start/end coordinates and orientation.
       - plot_layout.tsv: one row per sequence (reference and query) with
         order, length, cumulative offset, and centre position for axis
         tick placement and gridlines.

The greedy reordering is the same heuristic used by visualisation tools
such as D-GENIES. It does not guarantee a global optimum, but in practice
it places dominant orthologue pairs on the diagonal cleanly.
"""
from pathlib import Path

import pandas as pd

coords_path = Path(snakemake.input.coords)
ref_summary_path = Path(snakemake.input.ref_summary)
query_summary_path = Path(snakemake.input.query_summary)
out_coords = Path(snakemake.output.table)
out_layout = Path(snakemake.output.layout)

min_len = float(snakemake.params.min_len)
min_identity = float(snakemake.params.min_identity)

out_coords.parent.mkdir(parents=True, exist_ok=True)

cols = [
    "S1", "E1", "S2", "E2",
    "LEN1", "LEN2",
    "IDY",
    "LENR", "LENQ",
    "COVR", "COVQ",
    "REF", "QRY",
]

df = pd.read_csv(coords_path, sep="\t", header=None, comment="#")
if df.empty:
    raise SystemExit(f"No alignments found in {coords_path}")
if df.shape[1] < 13:
    raise SystemExit(
        f"Expected at least 13 columns from show-coords -rclTH, "
        f"found {df.shape[1]}"
    )
df = df.iloc[:, :13]
df.columns = cols

num_cols = ["S1", "E1", "S2", "E2", "LEN1", "LEN2",
            "IDY", "LENR", "LENQ", "COVR", "COVQ"]
for col in num_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.dropna(subset=num_cols + ["REF", "QRY"]).copy()

ref_summary = pd.read_csv(ref_summary_path, sep="\t")
query_summary = pd.read_csv(query_summary_path, sep="\t")

ref_order = ref_summary["sequence_id"].tolist()
ref_len = dict(zip(ref_summary["sequence_id"], ref_summary["length"]))
query_len = dict(zip(query_summary["sequence_id"], query_summary["length"]))

df = df[df["REF"].isin(ref_order) & df["QRY"].isin(query_len)].copy()

df["alignment_len_min"] = df[["LEN1", "LEN2"]].min(axis=1)
df = df[
    (df["alignment_len_min"] >= min_len) &
    (df["IDY"] >= min_identity)
].copy()

if df.empty:
    raise SystemExit(
        f"No alignments left after filtering for minimum length "
        f">= {min_len} bp and identity >= {min_identity}%"
    )

agg = (
    df.groupby(["REF", "QRY"])["alignment_len_min"]
      .sum()
      .reset_index()
      .rename(columns={"alignment_len_min": "total_aligned"})
)

assigned = []
remaining = list(query_summary["sequence_id"])

for ref_seq in ref_order:
    candidates = agg[
        (agg["REF"] == ref_seq) & (agg["QRY"].isin(remaining))
    ]
    if candidates.empty:
        continue
    best = candidates.sort_values(
        "total_aligned", ascending=False
    ).iloc[0]["QRY"]
    assigned.append(best)
    remaining.remove(best)

query_order = assigned + remaining


def cumulative_offsets(seq_ids, length_lookup):
    offsets = {}
    cursor = 0
    for sid in seq_ids:
        offsets[sid] = cursor
        cursor += length_lookup[sid]
    return offsets, cursor


ref_offset, ref_total = cumulative_offsets(ref_order, ref_len)
query_offset, query_total = cumulative_offsets(query_order, query_len)

df["x_start"] = df["REF"].map(ref_offset) + df["S1"]
df["x_end"]   = df["REF"].map(ref_offset) + df["E1"]
df["y_start"] = df["QRY"].map(query_offset) + df["S2"]
df["y_end"]   = df["QRY"].map(query_offset) + df["E2"]

df["orientation"] = (df["S2"] <= df["E2"]).map(
    {True: "Forward", False: "Reverse"}
)

keep_cols = [
    "REF", "QRY",
    "S1", "E1", "S2", "E2",
    "LEN1", "LEN2", "IDY",
    "orientation",
    "x_start", "x_end", "y_start", "y_end",
]
df_out = df[keep_cols].sort_values(["QRY", "REF", "S1", "S2"])
df_out.to_csv(out_coords, sep="\t", index=False)

layout_rows = []
for i, sid in enumerate(ref_order, start=1):
    layout_rows.append({
        "side": "reference",
        "order": i,
        "sequence_id": sid,
        "length": ref_len[sid],
        "offset": ref_offset[sid],
        "center": ref_offset[sid] + ref_len[sid] / 2,
    })
for i, sid in enumerate(query_order, start=1):
    layout_rows.append({
        "side": "query",
        "order": i,
        "sequence_id": sid,
        "length": query_len[sid],
        "offset": query_offset[sid],
        "center": query_offset[sid] + query_len[sid] / 2,
    })
pd.DataFrame(layout_rows).to_csv(out_layout, sep="\t", index=False)
