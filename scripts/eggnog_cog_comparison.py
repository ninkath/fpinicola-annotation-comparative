#!/usr/bin/env python3

"""
eggnog_cog_comparison.py

Compare COG functional category profiles across species
from eggNOG-mapper output files.

Produces:
  - A TSV table with COG category counts per species
  - A text summary with totals and annotation rates
"""

import csv
from collections import Counter
from pathlib import Path

# COG category descriptions
COG_DESCRIPTIONS = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


def parse_eggnog(filepath):
    """
    Parse eggNOG-mapper annotations file.
    Returns (total_proteins, annotated_count, cog_counter).
    COG categories are in column index 6 (COG_category).
    Each protein can have multiple COG categories (e.g., 'KT').
    """
    total = 0
    annotated = 0
    cog_counts = Counter()

    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            cols = line.strip().split("\t")
            if len(cols) > 6 and cols[6] and cols[6] != "-":
                annotated += 1
                # Each character is a separate COG category
                for cat in cols[6]:
                    if cat in COG_DESCRIPTIONS:
                        cog_counts[cat] += 1

    return total, annotated, cog_counts


def main():
    # Input files from snakemake
    species_files = {
        "F_pinicola": snakemake.input.pin,
        "F_schrenkii": snakemake.input.sch,
        "F_rosea": snakemake.input.ros,
    }

    results = {}
    for sp, fpath in species_files.items():
        total, annotated, cog_counts = parse_eggnog(fpath)
        results[sp] = {
            "total": total,
            "annotated": annotated,
            "cog_counts": cog_counts,
        }

    # Collect all observed COG categories
    all_cats = sorted(
        set().union(*(r["cog_counts"].keys() for r in results.values()))
    )
    species_names = list(species_files.keys())

    # Write comparison table
    with open(snakemake.output.table, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["COG_category", "Description"] + species_names
        writer.writerow(header)
        for cat in all_cats:
            row = [
                cat,
                COG_DESCRIPTIONS.get(cat, "Unknown"),
            ] + [results[sp]["cog_counts"].get(cat, 0) for sp in species_names]
            writer.writerow(row)

        # Add totals row
        writer.writerow([])
        writer.writerow(
            ["TOTAL_PROTEINS", ""]
            + [results[sp]["total"] for sp in species_names]
        )
        writer.writerow(
            ["COG_ANNOTATED", ""]
            + [results[sp]["annotated"] for sp in species_names]
        )
        writer.writerow(
            ["ANNOTATION_RATE", ""]
            + [
                f"{results[sp]['annotated']/results[sp]['total']*100:.1f}%"
                if results[sp]["total"] > 0
                else "0%"
                for sp in species_names
            ]
        )

    # Write text summary
    with open(snakemake.output.summary, "w") as f:
        f.write("COG Category Comparison Summary\n")
        f.write("=" * 60 + "\n\n")

        for sp in species_names:
            r = results[sp]
            rate = r["annotated"] / r["total"] * 100 if r["total"] > 0 else 0
            f.write(f"{sp}:\n")
            f.write(f"  Total proteins: {r['total']}\n")
            f.write(f"  With COG annotation: {r['annotated']} ({rate:.1f}%)\n")
            f.write(f"  Total COG assignments: {sum(r['cog_counts'].values())}\n")
            f.write("\n")

        f.write("\nCategories with largest between-species differences:\n")
        f.write("-" * 60 + "\n")

        # Find categories with largest variation
        diffs = []
        for cat in all_cats:
            vals = [results[sp]["cog_counts"].get(cat, 0) for sp in species_names]
            max_diff = max(vals) - min(vals)
            diffs.append((max_diff, cat, vals))

        diffs.sort(reverse=True)
        for diff, cat, vals in diffs[:10]:
            desc = COG_DESCRIPTIONS.get(cat, "Unknown")
            counts = ", ".join(
                f"{sp}: {v}" for sp, v in zip(species_names, vals)
            )
            f.write(f"  [{cat}] {desc}\n")
            f.write(f"      {counts} (max diff: {diff})\n")


if __name__ == "__main__":
    main()
