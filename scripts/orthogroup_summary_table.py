#!/usr/bin/env python3

from pathlib import Path
import csv

presence_path = Path(snakemake.input.presence)
single_copy_path = Path(snakemake.input.single_copy)
output_path = Path(snakemake.output.table)

output_path.parent.mkdir(parents=True, exist_ok=True)

species = ["F_pinicola", "F_schrenkii", "F_rosea"]

counts = {
    "total_orthogroups": 0,
    "shared_all_three": 0,
    "shared_F_pinicola_F_schrenkii_only": 0,
    "shared_F_pinicola_F_rosea_only": 0,
    "shared_F_schrenkii_F_rosea_only": 0,
    "unique_F_pinicola": 0,
    "unique_F_schrenkii": 0,
    "unique_F_rosea": 0,
}

with presence_path.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        counts["total_orthogroups"] += 1

        p = int(row["F_pinicola"])
        s = int(row["F_schrenkii"])
        r = int(row["F_rosea"])

        pattern = (p, s, r)

        if pattern == (1, 1, 1):
            counts["shared_all_three"] += 1
        elif pattern == (1, 1, 0):
            counts["shared_F_pinicola_F_schrenkii_only"] += 1
        elif pattern == (1, 0, 1):
            counts["shared_F_pinicola_F_rosea_only"] += 1
        elif pattern == (0, 1, 1):
            counts["shared_F_schrenkii_F_rosea_only"] += 1
        elif pattern == (1, 0, 0):
            counts["unique_F_pinicola"] += 1
        elif pattern == (0, 1, 0):
            counts["unique_F_schrenkii"] += 1
        elif pattern == (0, 0, 1):
            counts["unique_F_rosea"] += 1

single_copy_count = 0
with single_copy_path.open() as fh:
    for line in fh:
        line = line.strip()
        if line:
            single_copy_count += 1

counts["single_copy_orthologues"] = single_copy_count

ordered_metrics = [
    "total_orthogroups",
    "shared_all_three",
    "shared_F_pinicola_F_schrenkii_only",
    "shared_F_pinicola_F_rosea_only",
    "shared_F_schrenkii_F_rosea_only",
    "unique_F_pinicola",
    "unique_F_schrenkii",
    "unique_F_rosea",
    "single_copy_orthologues",
]

with output_path.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "count"])
    for metric in ordered_metrics:
        writer.writerow([metric, counts[metric]])
