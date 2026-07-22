#!/usr/bin/env python3
"""Reconcile the answer key's `consequence` field to the REAL SnpEff annotation.

WHY. An independent audit (run with the vcf2report ACMG engine) compared the manifest's `consequence`
column against what SnpEff (`GRCh38.mane.1.5.refseq`) actually calls at each planted coordinate in the
shipped `realistic_annotated/` VCFs. The coordinates, alleles, genes and zygosities were always correct
(0 mismatches) — but the `consequence` label was not: it carried a naive / source-derived value that
disagreed with the real annotator in 26 cases in `cohort.tsv` and 12 in `planted_variants.tsv`
(9 in both). Examples: "frameshift_variant" on the non-coding RNAs RNU5B-1 / RNU4-2 (no reading frame);
"missense_variant" for variants that are actually at a splice donor/acceptor site or intronic;
"stop_gained" for a genuine frameshift indel. A benchmark used to score an ANNOTATOR must carry the
consequence the real annotator produces, so this rewrites the field to ground truth.

WHAT. For the PRIMARY planted allele of every case, replace `consequence` with the most-severe SO term
SnpEff assigns to the planted gene at that coordinate, in: manifest/cohort.tsv,
manifest/planted_variants.tsv (primary rows only), and sidecars/SYN-NNN.planted.tsv (primary rows).
Only the `consequence` field changes; every other value is preserved byte-for-byte (line-based edit).
Idempotent: re-running it is a no-op once reconciled.

  python3 scripts/reconcile_consequence.py --annotated realistic_annotated [--check]
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _find_vcf(ann_dir: Path, sid: str) -> str | None:
    for p in glob.glob(f"{ann_dir}/**/{sid}.annotated.vcf.gz", recursive=True):
        if not os.path.basename(p).startswith("._"):  # skip macOS AppleDouble junk
            return p
    return None


def _most_severe_consequence(vcf: str, chrom: str, pos: str, ref: str, alt: str, gene: str) -> str | None:
    """The most-severe SO term SnpEff assigns to `gene` at this coordinate (SnpEff orders ANN by
    severity, so the first gene-matching entry's first &-term is the canonical consequence)."""
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom and f[1] == pos and f[3] == ref and alt in f[4].split(","):
                ann = [x[4:] for x in f[7].split(";") if x.startswith("ANN=")]
                if not ann:
                    return None
                for entry in ann[0].split(","):
                    p = entry.split("|")
                    if len(p) > 3 and p[3] == gene:
                        return p[1].split("&")[0]
                return ann[0].split(",")[0].split("|")[1].split("&")[0]  # gene fallback
    return None


def _primary_targets(ann_dir: Path) -> dict[str, str]:
    """syn_id -> real most-severe consequence, for the primary planted allele."""
    pv = REPO / "manifest" / "planted_variants.tsv"
    out: dict[str, str] = {}
    with open(pv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("allele") != "primary":
                continue
            vcf = _find_vcf(ann_dir, r["syn_id"])
            if not vcf:
                continue
            c = _most_severe_consequence(vcf, r["chrom"], r["pos"], r["ref"], r["alt"], r["gene"])
            if c:
                out[r["syn_id"]] = c
    return out


def _rewrite(path: Path, real: dict[str, str], primary_only: bool, check: bool) -> list[tuple]:
    """Line-based rewrite of the `consequence` column; returns [(syn_id, old, new), ...] changed."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    ci = header.index("consequence")
    ai = header.index("allele") if "allele" in header else None
    changed = []
    for i in range(1, len(lines)):
        cols = lines[i].split("\t")
        if len(cols) <= ci:
            continue
        sid = cols[0]
        if primary_only and ai is not None and (len(cols) <= ai or cols[ai] != "primary"):
            continue
        new = real.get(sid)
        if new and cols[ci] != new:
            changed.append((sid, cols[ci], new))
            cols[ci] = new
            lines[i] = "\t".join(cols)
    if changed and not check:
        path.write_text("\n".join(lines) + "\n")
    return changed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Reconcile the answer key consequence to real SnpEff.")
    ap.add_argument("--annotated", default="realistic_annotated", help="dir of SYN-NNN.annotated.vcf.gz")
    ap.add_argument("--check", action="store_true", help="report changes without writing")
    args = ap.parse_args(argv)

    ann_dir = Path(args.annotated)
    if not ann_dir.is_absolute():
        ann_dir = REPO / ann_dir
    real = _primary_targets(ann_dir)
    if not real:
        print(f"No annotated VCFs found under {ann_dir}. Run fetch.sh first.")
        return 1

    targets = [REPO / "manifest" / "cohort.tsv", REPO / "manifest" / "planted_variants.tsv"]
    targets += sorted((REPO / "sidecars").glob("SYN-*.planted.tsv"))
    total = 0
    for t in targets:
        primary_only = t.name != "cohort.tsv"
        ch = _rewrite(t, real, primary_only, args.check)
        total += len(ch)
        if ch:
            rel = t.relative_to(REPO)
            print(f"{'would change' if args.check else 'changed'} {len(ch):>2} in {rel}")
            for sid, old, new in ch[:6]:
                print(f"    {sid}: {old} -> {new}")
            if len(ch) > 6:
                print(f"    ... (+{len(ch)-6} more)")
    print(f"\n{'(check) ' if args.check else ''}total field changes: {total} across {len(targets)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
