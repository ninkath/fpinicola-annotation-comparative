#!/usr/bin/env python3

from pathlib import Path
import re

def read_fasta_stats(fasta_path: str):
    lengths = []
    gc_count = 0
    total_bases = 0

    seq_len = 0
    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                    seq_len = 0
            else:
                seq = line.upper()
                seq_len += len(seq)
                total_bases += len(seq)
                gc_count += seq.count("G") + seq.count("C")
        if seq_len > 0:
            lengths.append(seq_len)

    lengths_sorted = sorted(lengths, reverse=True)
    total_len = sum(lengths_sorted)

    half = total_len / 2
    running = 0
    n50 = 0
    for length in lengths_sorted:
        running += length
        if running >= half:
            n50 = length
            break

    gc_percent = round((gc_count / total_bases) * 100, 2) if total_bases > 0 else 0.0

    return {
        "assembly_size_bp": total_len,
        "num_sequences": len(lengths_sorted),
        "gc_percent": gc_percent,
        "N50_bp": n50,
    }


def read_busco_summary(busco_path: str):
    summary_line = None
    n_buscos = None

    with open(busco_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("C:"):
                summary_line = line
            elif line.endswith("Total BUSCO groups searched"):
                parts = line.split()
                try:
                    n_buscos = int(parts[0])
                except Exception:
                    pass

    if summary_line is None:
        raise ValueError(f"Could not find BUSCO summary line in {busco_path}")

    m = re.search(
        r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)",
        summary_line
    )
    if not m:
        raise ValueError(f"Could not parse BUSCO summary line in {busco_path}: {summary_line}")

    return {
        "busco_complete_pct": float(m.group("C")),
        "busco_single_pct": float(m.group("S")),
        "busco_duplicated_pct": float(m.group("D")),
        "busco_fragmented_pct": float(m.group("F")),
        "busco_missing_pct": float(m.group("M")),
        "busco_n": int(m.group("n")) if m.group("n") else n_buscos,
    }


species_info = [
    ("F_pinicola", snakemake.input.pin_assembly, snakemake.input.pin_busco),
    ("F_schrenkii", snakemake.input.sch_assembly, snakemake.input.sch_busco),
    ("F_rosea", snakemake.input.ros_assembly, snakemake.input.ros_busco),
]

out_path = Path(snakemake.output.table)
out_path.parent.mkdir(parents=True, exist_ok=True)

rows = []
for species, assembly_path, busco_path in species_info:
    fasta_stats = read_fasta_stats(assembly_path)
    busco_stats = read_busco_summary(busco_path)

    row = {
        "species": species,
        **fasta_stats,
        **busco_stats,
    }
    rows.append(row)

header = [
    "species",
    "assembly_size_bp",
    "num_sequences",
    "N50_bp",
    "gc_percent",
    "busco_complete_pct",
    "busco_single_pct",
    "busco_duplicated_pct",
    "busco_fragmented_pct",
    "busco_missing_pct",
    "busco_n",
]

with open(out_path, "w") as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(str(row[col]) for col in header) + "\n")
