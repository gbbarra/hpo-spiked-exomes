#!/usr/bin/env bash
# Annotate a raw exome VCF with gene + consequence + HGVS (SnpEff) — the benchmark's annotator.
#
# The answer key's `consequence` is defined as whatever THIS step produces (read back by
# scripts/reconcile_consequence.py). A raw DRAGEN VCF carries only quality fields; SnpEff adds
# INFO/ANN with the molecular consequence on MANE transcripts.
#
#   bash scripts/build/annotate.sh RAW.vcf.gz OUT.annotated.vcf.gz [REF_GRCh38.fa]
#
# SnpEff is resolved from $SNPEFF_JAR, else `snpEff` on PATH. DB defaults to GRCh38.mane.1.5.refseq
# (override with $SNPEFF_DB) — this MUST match the DB the shipped realistic_annotated/ was built with,
# or the reconciled consequences won't be comparable across releases.
set -euo pipefail

RAW="${1:?raw VCF (bgzipped) required}"
OUT="${2:?output path required}"
REF="${3:-}"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.mane.1.5.refseq}"
THREADS="${THREADS:-4}"

if [[ -n "${SNPEFF_JAR:-}" ]]; then SNPEFF=(java -Xmx8g -jar "$SNPEFF_JAR")
elif command -v snpEff >/dev/null 2>&1; then SNPEFF=(snpEff)
else echo "ERROR: SnpEff not found (set \$SNPEFF_JAR or put snpEff on PATH)." >&2; exit 1; fi
for t in bcftools bgzip tabix; do command -v "$t" >/dev/null || { echo "ERROR: '$t' not found." >&2; exit 1; }; done

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# Chromosome naming must match the DB or SnpEff silently annotates NOTHING (every record becomes
# ERROR_CHROMOSOME_NOT_FOUND). MANE GRCh38 DBs are Ensembl-style ("1"); DRAGEN is UCSC-style ("chr1").
# Rename in, and back out, so $OUT keeps the caller's naming.
vcf_chr="$(set +o pipefail; bcftools view -H "$RAW" 2>/dev/null | head -1 | cut -f1)"
db_chr="$("${SNPEFF[@]}" genes2bed "$SNPEFF_DB" TP53 2>/dev/null | sed -n '2p' | cut -f1)"
vcf_ucsc=0; [[ "$vcf_chr" == chr* ]] && vcf_ucsc=1
db_ucsc=0;  [[ "$db_chr"  == chr* ]] && db_ucsc=1
[[ -z "$db_chr" ]] && db_ucsc="${SNPEFF_DB_UCSC:-0}"
rename_needed=0; [[ "$vcf_ucsc" != "$db_ucsc" ]] && rename_needed=1
if [[ "$rename_needed" == 1 ]]; then
  for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > "$tmp/to_ens.txt"; echo "chrM MT" >> "$tmp/to_ens.txt"
  for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > "$tmp/to_ucsc.txt"; echo "MT chrM" >> "$tmp/to_ucsc.txt"
  if [[ "$vcf_ucsc" == 1 ]]; then IN_MAP="$tmp/to_ens.txt"; OUT_MAP="$tmp/to_ucsc.txt"
  else IN_MAP="$tmp/to_ucsc.txt"; OUT_MAP="$tmp/to_ens.txt"; fi
fi

echo "[1/2] normalizing (split multiallelics${REF:+ + left-align}) ..." >&2
norm_args=(-m -any --threads "$THREADS"); [[ -n "$REF" ]] && norm_args+=(-f "$REF")
bcftools norm "${norm_args[@]}" -Oz -o "$tmp/norm.vcf.gz" "$RAW"

echo "[2/2] SnpEff consequence + HGVS ($SNPEFF_DB) ..." >&2
if [[ "$rename_needed" == 1 ]]; then
  bcftools annotate --rename-chrs "$IN_MAP" "$tmp/norm.vcf.gz" -Oz -o "$tmp/in.vcf.gz"
else mv "$tmp/norm.vcf.gz" "$tmp/in.vcf.gz"; fi
# No -canon: with a MANE DB it would drop MANE Plus Clinical (a false negative for exactly the
# genes that need it). SnpEff orders ANN by severity, so ANN[0] stays the most-severe consequence.
"${SNPEFF[@]}" -noStats -hgvs "$SNPEFF_DB" "$tmp/in.vcf.gz" > "$tmp/snpeff.vcf"
if [[ "$rename_needed" == 1 ]]; then
  bcftools annotate --rename-chrs "$OUT_MAP" "$tmp/snpeff.vcf" -Oz -o "$tmp/ann.vcf.gz"
else bgzip -c "$tmp/snpeff.vcf" > "$tmp/ann.vcf.gz"; fi

# Fail loudly on the silent-failure mode (naming/DB mismatch -> a VCF that looks fine but has no ANN).
total="$(bcftools view -H "$tmp/ann.vcf.gz" | wc -l | tr -d ' ')"
with_ann="$(bcftools view -H "$tmp/ann.vcf.gz" | grep -c 'ANN=' || true)"
if [[ "$total" -gt 0 && "$with_ann" -lt $((total / 2)) ]]; then
  echo "ERROR: only $with_ann/$total records got ANN — annotation failed (VCF=$vcf_chr db=$db_chr)." >&2
  exit 1
fi
echo "      annotated $with_ann/$total records" >&2

bcftools index -t "$tmp/ann.vcf.gz"
mkdir -p "$(dirname "$OUT")"
mv -f "$tmp/ann.vcf.gz" "$OUT"; mv -f "$tmp/ann.vcf.gz.tbi" "$OUT.tbi"
echo "Done -> $OUT" >&2
