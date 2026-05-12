#!/usr/bin/env python3
"""
Summary statistics for pairwise NUCmer dotplots.
 
Reads the per-comparison plot_coords.tsv and plot_layout.tsv files and
produces a single summary table covering all three pairwise alignments.
 
Run from the project root after the synteny dotplots have been generated:
 
    python scripts/synteny_summary_table.py
"""
from pathlib import Path
import pandas as pd
 
base = Path("results/comparative/synteny_dotplot")
 
comparisons = [
    ("F_rosea_vs_F_pinicola",     "F. rosea vs F. pinicola"),
    ("F_schrenkii_vs_F_pinicola", "F. schrenkii vs F. pinicola"),
    ("F_rosea_vs_F_schrenkii",    "F. rosea vs F. schrenkii"),
]
 
rows = []
for slug, label in comparisons:
    coords_path = base / slug / "plot_coords.tsv"
    layout_path = base / slug / "plot_layout.tsv"
    coords = pd.read_csv(coords_path, sep="\t")
    layout = pd.read_csv(layout_path, sep="\t")
    ref_layout = layout[layout["side"] == "reference"]
    qry_layout = layout[layout["side"] == "query"]
    ref_len = ref_layout["length"].sum()
    rows.append({
        "Comparison": label,
        "Ref sequences": len(ref_layout),
        "Query sequences": len(qry_layout),
        "Retained alignments": len(coords),
        "Mean identity (%)": round(coords["IDY"].mean(), 2),
        "Total aligned on ref (Mb)": round(coords["LEN1"].sum() / 1e6, 2),
        "Ref covered (%)": round(coords["LEN1"].sum() / ref_len * 100, 2),
    })
 
df = pd.DataFrame(rows)
out = base / "alignment_summary.tsv"
df.to_csv(out, sep="\t", index=False)
 
print(df.to_string(index=False))
print(f"\nWrote {out}")
 
