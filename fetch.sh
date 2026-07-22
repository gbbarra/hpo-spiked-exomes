#!/usr/bin/env bash
# Download the hpo-spiked-exomes VCFs (release assets), verify checksums, and extract them into
# ./realistic (raw tell-free) and ./realistic_annotated (SnpEff-annotated).
#   bash fetch.sh [TAG]        # TAG defaults to data-v1
set -euo pipefail
REPO="gbbarra/hpo-spiked-exomes"
TAG="${1:-data-v1}"
cd "$(cd "$(dirname "$0")" && pwd)"

command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }

echo "Fetching $REPO release $TAG ..." >&2
if command -v gh >/dev/null; then
  gh release download "$TAG" -R "$REPO" -p '*.tar.zst' -p 'SHA256SUMS' --clobber
else
  base="https://github.com/$REPO/releases/download/$TAG"
  for a in hpo_spiked_exomes_realistic.tar.zst hpo_spiked_exomes_realistic_annotated.tar.zst SHA256SUMS; do
    curl -fL --retry 3 -o "$a" "$base/$a"
  done
fi

echo "Verifying checksums ..." >&2
shasum -a 256 -c SHA256SUMS

mkdir -p realistic realistic_annotated
zstd -dc hpo_spiked_exomes_realistic.tar.zst           | tar -x -C realistic
zstd -dc hpo_spiked_exomes_realistic_annotated.tar.zst | tar -x -C realistic_annotated
echo "Done: $(ls realistic/SYN-*.vcf.gz 2>/dev/null | wc -l | tr -d ' ') raw + $(ls realistic_annotated/SYN-*.annotated.vcf.gz 2>/dev/null | wc -l | tr -d ' ') annotated exomes." >&2
