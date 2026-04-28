#!/usr/bin/env python3
"""
parse_repeatmasker_out.py

Parse RepeatMasker .out file into a BED with simplified repeat classes
for circos visualization.

Classes:
  LTR          - All LTR retrotransposons (Ty1/Copia, Gypsy/DIRS1, etc.)
  DNA_TE       - All DNA transposons (hobo-Activator, Tc1, rolling-circles, etc.)
  LINE         - LINEs (L1/CIN4, etc.)
  Unclassified - Unknown/unclassified interspersed repeats
  Other        - Small RNA, simple repeats, low complexity, satellites, etc.

Output BED columns: chrom, start, end, class, repeat_name, divergence

Usage:
    python3 parse_repeatmasker_out.py \
        --input assembly.fasta.out \
        --output repeats_by_class.bed
"""

import argparse
import sys


def classify_repeat(class_family):
    """
    Map RepeatMasker class/family string to simplified class.
    """
    cf = class_family.strip()
    cfl = cf.lower()

    # LTR retrotransposons
    if "ltr" in cfl:
        return "LTR"

    # DNA transposons
    if cfl.startswith("dna") or "dna/" in cfl:
        return "DNA_TE"
    if any(k in cfl for k in ["hobo-activator", "tc1-is630-pogo",
                               "piggybac", "tourist", "helitron",
                               "rolling-circle", "rc/"]):
        return "DNA_TE"

    # LINEs
    if cfl.startswith("line") or "line/" in cfl:
        return "LINE"

    # SINEs (rare in fungi but handle anyway)
    if cfl.startswith("sine") or "sine/" in cfl:
        return "Other"

    # Other retroelements not caught above
    if "retroposon" in cfl:
        return "LINE"

    # Unclassified interspersed
    if "unknown" in cfl or "unspecified" in cfl or cf == "Unknown":
        return "Unclassified"

    # Simple repeats, low complexity, satellites, rRNA, snRNA, etc.
    if any(k in cfl for k in ["simple_repeat", "low_complexity",
                               "satellite", "rrna", "snrna",
                               "scrna", "trna", "srpna",
                               "small_rna", "other"]):
        return "Other"

    # Catch-all for anything not classified above
    return "Unclassified"


def parse_rm_out(input_path, output_path):
    """
    Parse RepeatMasker .out file.
    The .out format has 3 header lines, then data columns:
      0: SW score
      1: %divergence
      2: %deletions
      3: %insertions
      4: query sequence (contig)
      5: query begin
      6: query end
      7: query (left)
      8: strand (+ or C)
      9: repeat name
     10: class/family
     11-14: repeat position info
    """
    n_parsed = 0
    n_skipped = 0
    class_counts = {}

    with open(input_path) as fin, open(output_path, "w") as fout:
        for i, line in enumerate(fin):
            # Skip header lines (first 3 lines)
            if i < 3:
                continue

            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 11:
                n_skipped += 1
                continue

            try:
                divergence = float(parts[1])
                contig = parts[4]
                start = int(parts[5]) - 1  # convert to 0-based BED
                end = int(parts[6])
                repeat_name = parts[9]
                class_family = parts[10]
            except (ValueError, IndexError):
                n_skipped += 1
                continue

            simplified = classify_repeat(class_family)

            fout.write(f"{contig}\t{start}\t{end}\t{simplified}\t"
                       f"{repeat_name}\t{divergence:.1f}\n")

            class_counts[simplified] = class_counts.get(simplified, 0) + 1
            n_parsed += 1

    # Summary to stderr
    print(f"Parsed {n_parsed} repeat elements, skipped {n_skipped} lines",
          file=sys.stderr)
    for cls in sorted(class_counts.keys()):
        print(f"  {cls}: {class_counts[cls]}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True,
                        help="RepeatMasker .out file")
    parser.add_argument("--output", required=True,
                        help="Output BED file with repeat classes")
    args = parser.parse_args()

    parse_rm_out(args.input, args.output)
