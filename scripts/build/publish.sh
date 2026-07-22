#!/usr/bin/env bash
# Pack the tell-free realistic + annotated VCFs into the release tarballs and (optionally) upload.
#
#   bash scripts/build/publish.sh <realistic_dir> <annotated_dir> <TAG> [--upload]
#
# Produces, under a work dir:
#   hpo_spiked_exomes_realistic.tar.zst            (raw, tell-free)
#   hpo_spiked_exomes_realistic_annotated.tar.zst  (SnpEff-annotated)
#   SHA256SUMS                                     (shasum -a 256 — the format fetch.sh verifies)
#
# COPYFILE_DISABLE=1 is the whole point: without it, macOS bsdtar synthesizes `._*` AppleDouble
# sidecars from xattrs (a `--exclude='._*'` does NOT stop that — the files aren't on disk).
set -euo pipefail
export COPYFILE_DISABLE=1

REAL="${1:?realistic dir required}"
ANN="${2:?annotated dir required}"
TAG="${3:?release tag required (e.g. data-v2)}"
UPLOAD=0; [ "${4:-}" = "--upload" ] && UPLOAD=1
REPO_SLUG="gbbarra/hpo-spiked-exomes"

command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }
[ "$UPLOAD" = 1 ] && { command -v gh >/dev/null || { echo "ERROR: gh not found." >&2; exit 1; }; }
[ -d "$REAL" ] && [ -d "$ANN" ] || { echo "ERROR: dirs not found." >&2; exit 1; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
R="hpo_spiked_exomes_realistic.tar.zst"
A="hpo_spiked_exomes_realistic_annotated.tar.zst"

echo ">>> packing raw tell-free ($(find "$REAL" -name 'SYN-*.vcf.gz' | wc -l | tr -d ' ') VCFs) ..." >&2
( cd "$REAL" && tar --exclude='._*' --exclude='.DS_Store' -cf - SYN-*.vcf.gz ) | zstd -19 -T0 -q -o "$WORK/$R"
echo ">>> packing annotated ($(find "$ANN" -name 'SYN-*.annotated.vcf.gz' | wc -l | tr -d ' ') VCFs) ..." >&2
( cd "$ANN" && tar --exclude='._*' --exclude='.DS_Store' -cf - SYN-*.annotated.vcf.gz SYN-*.annotated.vcf.gz.tbi ) \
    | zstd -19 -T0 -q -o "$WORK/$A"

# sanity: 0 ._* in either tarball
for f in "$R" "$A"; do
  n="$(zstd -dc "$WORK/$f" | tar -tf - | grep -cE '(^|/)\._' || true)"
  [ "$n" = 0 ] || { echo "ERROR: $n ._* entries in $f" >&2; exit 1; }
done
( cd "$WORK" && shasum -a 256 "$R" "$A" > SHA256SUMS && sed 's/^/    /' SHA256SUMS >&2 )

if [ "$UPLOAD" = 1 ]; then
  echo ">>> uploading to $TAG ..." >&2
  ( cd "$WORK" && gh release upload "$TAG" -R "$REPO_SLUG" --clobber "$R" "$A" SHA256SUMS )
  echo ">>> uploaded." >&2
else
  echo ">>> built (NOT uploaded — pass --upload). Artifacts in $WORK" >&2
  echo "$WORK"
fi
