#!/usr/bin/env bash
# Build the raw spiked VCFs for a cohort TSV: for each row, stream that 1000G DRAGEN v4.4.7 exome,
# subset to the MANE/GENCODE exome BED, and plant the EXACT phenopacket variant. This plants a
# MARKER-BEARING record (SPIKED/GENE/CSQ/CLN*) on purpose — build_biallelic.py and
# realisticize_cohort.py read those markers, and realisticize is what strips them and captures the
# truth into the manifest/sidecars, producing the final tell-free VCF. Writes
# <out>/SYN-NNN.synthetic.vcf.gz + <out>/SYN-NNN.hpo.txt.
#
# SELF-CONTAINED: curl, bcftools, bgzip, tabix, python3. DRAGEN VCFs stream from the public S3 bucket.
#
#   BED=/path/exome_hg38.bed CLINVAR_VCF=/path/clinvar.vcf.gz \
#     bash scripts/build/make_cohort.sh new_cohort.tsv out_dir
#
# The MANE exome BED and the ClinVar GRCh38 VCF are NOT in this repo — point $BED / $CLINVAR_VCF at
# your copies (both are in the vcf2report build tree). Resumable: an existing SYN-NNN is skipped.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"; SCRIPTS="$(cd "$HERE/.." && pwd)"
COHORT="${1:?cohort TSV required}"
OUT="${2:?output dir required}"
BED="${BED:?set BED=/path/to/exome_hg38.bed}"
CLINVAR="${CLINVAR_VCF:?set CLINVAR_VCF=/path/to/clinvar_GRCh38.vcf.gz}"
BUCKET="1000genomes-dragen-v4-4-7"
PIPE="data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"
LIMIT="${N:-0}"

for t in curl bcftools bgzip tabix python3; do command -v "$t" >/dev/null || { echo "ERROR: '$t' not found." >&2; exit 1; }; done
for f in "$BED" "$CLINVAR" "$COHORT"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
mkdir -p "$OUT"

count=0
tail -n +2 "$COHORT" | while IFS=$'\t' read -r syn sample gene chrom pos ref alt cons disease hpo; do
  [ -n "${syn:-}" ] || continue
  count=$((count+1)); [ "$LIMIT" -gt 0 ] && [ "$count" -gt "$LIMIT" ] && break
  [ -f "$OUT/$syn.synthetic.vcf.gz" ] && { echo "==== [$syn] present — skip"; continue; }
  echo "==== [$syn] $sample $gene $chrom:$pos $ref>$alt ===="
  work="$OUT/.work/$syn"; mkdir -p "$work"

  url="https://$BUCKET.s3.amazonaws.com/$PIPE/$sample/$sample.hard-filtered.vcf.gz"
  if ! curl -fsI --max-time 30 "$url" >/dev/null 2>&1; then
    key="$(curl -s --max-time 30 "https://$BUCKET.s3.amazonaws.com/?list-type=2&prefix=$PIPE/$sample/&max-keys=200" \
           | grep -oE "<Key>[^<]+hard-filtered\.vcf\.gz</Key>" | sed 's/<[^>]*>//g' | head -1 || true)"
    [ -n "$key" ] || { echo "  WARN: no DRAGEN VCF for $sample — skip"; rm -rf "$work"; continue; }
    url="https://$BUCKET.s3.amazonaws.com/$key"
  fi
  if ! curl -fL --retry 6 --retry-delay 2 --speed-limit 1000000 --speed-time 25 \
            --connect-timeout 20 --max-time 3600 -C - -o "$work/raw.vcf.gz" "$url"; then
    echo "  WARN: $sample download failed — skip (resumable)"; rm -rf "$work"; continue
  fi

  # split multiallelics + subset to the MANE exome BED (one streaming pass)
  bcftools norm -m -any "$work/raw.vcf.gz" -Ou 2>/dev/null | bcftools view -T "$BED" -Oz -o "$work/exome.vcf.gz"

  # plant the exact phenopacket variant (marker-bearing; realisticize makes it tell-free later)
  python3 "$SCRIPTS/spike_variant.py" \
    --exome "$work/exome.vcf.gz" --clinvar "$CLINVAR" \
    --chrom "$chrom" --pos "$pos" --ref "$ref" --alt "$alt" --gene "$gene" \
    --consequence "$cons" --disease "$disease" --sample-id "$syn" --out "$OUT/$syn.synthetic.vcf"
  bgzip -f "$OUT/$syn.synthetic.vcf"; tabix -f -p vcf "$OUT/$syn.synthetic.vcf.gz"

  printf '%s\n' "$hpo" | tr ',' '\n' | sed '/^$/d' > "$OUT/$syn.hpo.txt"
  echo "  -> $OUT/$syn.synthetic.vcf.gz"
  rm -rf "$work"
done
rmdir "$OUT/.work" 2>/dev/null || true
echo "Done -> $OUT"
