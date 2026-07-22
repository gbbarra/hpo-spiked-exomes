#!/usr/bin/env python3
"""Apply the FAITHFUL genotype from the phenopacket plan to the raw spiked VCFs.

A lone het in an autosomal-recessive gene is a *carrier*, not a diagnosis. select_cases.py records
each case's real zygosity in the plan; this materialises it on the raw spiked VCF:

  * hom          — flip the planted call's genotype to 1/1.
  * compound_het — add the SECOND allele from the source case as a second heterozygous call.
  * single_het   — unchanged (the source genuinely recorded one allele).

The second allele carries only GENE + a SPIKED2 flag — NEVER a hand-written consequence, so it can
never disagree with the real SnpEff annotation (which is added on the next, annotate step). Truth for
the second allele is tracked externally by realisticize_cohort.py (the `second` rows in the manifest).

  python3 scripts/build/build_biallelic.py --plan new_plan.json --cohort-tsv new_cohort.tsv \
      --src-dir <raw vcf dir> --out <out dir>
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # for spike_pathogenic._CHROM_ORDER
from spike_pathogenic import _CHROM_ORDER  # noqa: E402


def _spike2_line(chrom, pos, ref, alt, gene, ncols):
    info = f"GENE={gene};SPIKED2=1"
    row = [chrom, str(pos), ".", ref, alt, "800", "PASS", info, "GT:DP:GQ:AD", "0/1:40:99:20,20"]
    return row + ["." for _ in range(ncols - len(row))]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", required=True)
    ap.add_argument("--cohort-tsv", required=True)
    ap.add_argument("--src-dir", required=True, help="dir with the raw <syn>.synthetic.vcf.gz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    plan = json.load(open(a.plan))
    rows = list(csv.DictReader(open(a.cohort_tsv), delimiter="\t"))
    src, out = Path(a.src_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    made = {"compound_het": 0, "hom": 0, "single_het": 0}

    for r in rows:
        syn, gene, pos = r["syn_id"], r["gene"], r["pos"]
        p = plan.get(syn) or {"mode": "single_het"}
        with gzip.open(src / f"{syn}.synthetic.vcf.gz", "rt") as fh:
            lines = fh.read().splitlines()
        meta = [l for l in lines if l.startswith("##")]
        header = next(l for l in lines if l.startswith("#CHROM"))
        ncols = len(header.split("\t"))
        body = [l.split("\t") for l in lines if l and not l.startswith("#")]

        if p["mode"] == "hom":
            for f in body:
                if f[1] == pos and "SPIKED=1" in f[7]:
                    g = f[9].split(":"); g[0] = "1/1"
                    if len(g) >= 4:
                        g[3] = "2,40"
                    f[9] = ":".join(g)
            made["hom"] += 1
        elif p["mode"] == "compound_het":
            body.append(_spike2_line(p["chrom"], p["pos2"], p["ref2"], p["alt2"], gene, ncols))
            body.sort(key=lambda f: (_CHROM_ORDER.get(f[0].replace("chr", ""), 99), int(f[1])))
            made["compound_het"] += 1
        else:
            made["single_het"] += 1

        extra = ['##INFO=<ID=SPIKED2,Number=0,Type=Flag,Description="Second biallelic spiked variant">']
        meta_out = meta + [m for m in extra if m not in meta]
        tmp = out / f"{syn}.synthetic.vcf"
        with open(tmp, "w") as w:
            w.write("\n".join(meta_out) + "\n" + header + "\n")
            for f in body:
                w.write("\t".join(f) + "\n")
        subprocess.run(["bgzip", "-f", str(tmp)], check=True)

    print(f"compound_het: {made['compound_het']} | hom: {made['hom']} | single_het: {made['single_het']}")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
