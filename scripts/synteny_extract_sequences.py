#!/usr/bin/env python3
"""
Extract assembled sequences for pairwise synteny analysis.

Sequences are filtered by minimum length, sorted in length-descending order,
and capped at a maximum count to keep downstream plots readable. The output
FASTA feeds into pairwise NUCmer alignment.

min_length can be either an integer (applied to all species) or a mapping
of species name to integer with an optional "default" key.
"""
from pathlib import Path
from Bio import SeqIO

input_fasta = Path(snakemake.input.fasta)
output_fasta = Path(snakemake.output.fasta)
summary_out = Path(snakemake.output.summary)
ml_param = snakemake.params.min_length
max_sequences = snakemake.params.max_sequences

# Resolve per-species override if a dict was passed.
if isinstance(ml_param, dict):
    species = snakemake.wildcards.species
    ml_value = ml_param.get(species, ml_param.get("default", 100_000))
else:
    ml_value = ml_param
min_length = int(ml_value)

output_fasta.parent.mkdir(parents=True, exist_ok=True)
summary_out.parent.mkdir(parents=True, exist_ok=True)

records = list(SeqIO.parse(input_fasta, "fasta"))
if not records:
    raise SystemExit(f"No sequences found in {input_fasta}")

records_sorted = sorted(records, key=lambda r: len(r.seq), reverse=True)
selected = [r for r in records_sorted if len(r.seq) >= min_length]

if not selected:
    raise SystemExit(
        f"No sequences in {input_fasta} reach the minimum length of "
        f"{min_length} bp"
    )

if max_sequences is not None and int(max_sequences) > 0:
    selected = selected[:int(max_sequences)]

SeqIO.write(selected, output_fasta, "fasta")

with open(summary_out, "w") as f:
    f.write("rank\tsequence_id\tlength\n")
    for i, rec in enumerate(selected, start=1):
        f.write(f"{i}\t{rec.id}\t{len(rec.seq)}\n")
