#!/usr/bin/env bash
set -euo pipefail

mkdir -p proteins_raw proteins_evidence

# ── 0. Download Swiss-Prot fungi subset (if not already present) ──
if [ ! -f proteins_raw/uniprot_sprot_fungi.fasta ]; then
  echo "Downloading Swiss-Prot fungi subset..."
  wget -q -O proteins_raw/uniprot_sprot_fungi.dat.gz \
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_fungi.dat.gz

  echo "Converting .dat to FASTA..."
  gunzip -k proteins_raw/uniprot_sprot_fungi.dat.gz
  awk '
    /^AC   / && ac=="" { ac=$2; gsub(/;/,"",ac) }
    /^DE   RecName: Full=/ && de=="" {
      de=$0; sub(/^DE   RecName: Full=/,"",de); sub(/ \{.*$/,"",de); sub(/;$/,"",de)
    }
    /^OS   / && os=="" { os=$0; sub(/^OS   /,"",os); sub(/\.$/,"",os) }
    /^SQ   / { reading=1; seq=""; next }
    reading && /^\/\// {
      printf ">sp|%s %s OS=%s\n%s\n", ac, de, os, seq
      ac=""; de=""; os=""; seq=""; reading=0
    }
    reading { gsub(/ /,"",$0); seq=seq $0 }
  ' proteins_raw/uniprot_sprot_fungi.dat > proteins_raw/uniprot_sprot_fungi.fasta
  rm -f proteins_raw/uniprot_sprot_fungi.dat
  echo "  Sequences: $(grep -c '^>' proteins_raw/uniprot_sprot_fungi.fasta)"
fi

# ── Helper: clean and assign safe IDs ──
# Usage: safe_rename PREFIX input.fa output.fa mapping.tsv
safe_rename () {
  local prefix="$1"
  local infile="$2"
  local outfile="$3"
  local mapfile="$4"
  echo "Processing: $infile -> $outfile (prefix: ${prefix})"
  awk -v prefix="$prefix" '
    BEGIN { n=0 }
    /^>/ {
      n++
      safeid = sprintf("%s%06d", prefix, n)
      # Store original header (without ">")
      orig = substr($0, 2)
      printf "%s\t%s\n", safeid, orig > "/dev/stderr"
      print ">" safeid
      next
    }
    { gsub(/\*$/, "", $0); print }
  ' "$infile" > "$outfile" 2>> "$mapfile"
}

# ── 1. Clean raw files with safe IDs ──
MAPFILE="proteins_evidence/id_mapping.tsv"
echo -e "safe_id\toriginal_header" > "$MAPFILE"

safe_rename "FBET" proteins_raw/f_betulina.faa \
  proteins_evidence/f_betulina.safe.fa "$MAPFILE"

safe_rename "FROS" proteins_raw/Fomitopsis_rosea_gca_004679265.ASM467926v1.pep.all.fa \
  proteins_evidence/f_rosea.safe.fa "$MAPFILE"

safe_rename "PPLA" proteins_raw/Postia_placenta_mad_698_r_gca_000006255.Postia_placenta_V1.0.pep.all.fa \
  proteins_evidence/p_placenta.safe.fa "$MAPFILE"

safe_rename "FPI3" proteins_raw/Fompi3_GeneCatalog_proteins_20120705.aa.fasta \
  proteins_evidence/fompi3.safe.fa "$MAPFILE"

safe_rename "SPFU" proteins_raw/uniprot_sprot_fungi.fasta \
  proteins_evidence/swissprot_fungi.safe.fa "$MAPFILE"

# ── 2. Create SET A (close relatives only) ──
cat \
  proteins_evidence/f_betulina.safe.fa \
  proteins_evidence/f_rosea.safe.fa \
  proteins_evidence/p_placenta.safe.fa \
  proteins_evidence/fompi3.safe.fa \
  > proteins_evidence/protein_evidence_SET_A_close_relatives.fa

# ── 3. Create SET B (SET A + Swiss-Prot fungi) ──
cat \
  proteins_evidence/protein_evidence_SET_A_close_relatives.fa \
  proteins_evidence/swissprot_fungi.safe.fa \
  > proteins_evidence/protein_evidence_SET_B_relatives_plus_swissprot_fungi.fa

# ── 4. Verify ──
echo ""
echo "Counts:"
echo -n "  SET_A: "; grep -c "^>" proteins_evidence/protein_evidence_SET_A_close_relatives.fa
echo -n "  SET_B: "; grep -c "^>" proteins_evidence/protein_evidence_SET_B_relatives_plus_swissprot_fungi.fa
echo -n "  Swiss-Prot fungi: "; grep -c "^>" proteins_evidence/swissprot_fungi.safe.fa
echo -n "  Mapping entries: "; tail -n +2 "$MAPFILE" | wc -l

echo ""
echo "Sample headers (SET_B, first 5):"
grep "^>" proteins_evidence/protein_evidence_SET_B_relatives_plus_swissprot_fungi.fa | head -5

echo ""
echo "Sample mapping (first 5):"
head -6 "$MAPFILE"

echo ""
echo "Checking for unsafe characters in headers [should be empty]:"
grep "^>" proteins_evidence/protein_evidence_SET_B_relatives_plus_swissprot_fungi.fa \
  | grep -E '[^>A-Za-z0-9]' -m 5 || echo "  (all clean)"

echo ""
echo "Invalid AA char check (SET_B) [should be empty]:"
grep -v "^>" proteins_evidence/protein_evidence_SET_B_relatives_plus_swissprot_fungi.fa \
  | grep -n "[^ACDEFGHIKLMNPQRSTVWYXBZUOJ]" -m 5 || true
