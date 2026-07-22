#!/usr/bin/env python3
"""Select NEW planted-case candidates from a GA4GH Phenopacket Store, for expanding the cohort.

Scans a phenopacket directory, picks cases whose causative variant + gene + coordinate are NOT
already in the shipped answer key, and emits the two build inputs:

  * a cohort TSV  (syn_id, sample, gene, chrom, pos, ref, alt, consequence, disease, hpo)
  * a faithful-genotype plan JSON  (per syn_id: hom / single_het / compound_het + 2nd allele)

New SYN ids continue after the existing maximum; new 1000G backgrounds are drawn from --backgrounds
and never reuse a sample already in the cohort.

CORRECTIONS baked in (see BUILD.md → "Pitfalls"):
  * disease is READ FROM THE PACKET and stored (never left blank).
  * the FULL, ordered HPO set is stored (never capped at 6; the sidecar must match it).
  * disease / HPO are written verbatim (never truncated — the 60-char bug).
  * the causative variant must be interpretationStatus == CAUSATIVE.
  * genes AND coordinates already in manifest/planted_variants.tsv are excluded (no reuse).
  * the `consequence` written here is only a COARSE PLACEHOLDER — the authoritative value is set
    later by scripts/reconcile_consequence.py from the real SnpEff annotation. Do NOT trust it.

  python3 scripts/build/select_cases.py --phenopackets /path/to/phenopacket-store \
      --backgrounds backgrounds.txt --n 100 --out-cohort new_cohort.tsv --out-plan new_plan.json
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenopacket import load_phenopacket  # noqa: E402

REPO = Path(__file__).resolve().parent.parent.parent


def _norm(c: str) -> str:
    return str(c).replace("chr", "").replace("CHR", "")


def _coarse_consequence(hgvs_p: str | None, ref: str, alt: str) -> str:
    """A PLACEHOLDER only — reconcile_consequence.py overwrites this with the real SnpEff term.
    Never used for scoring; it just seeds spike_variant's --consequence for coords absent from ClinVar."""
    hp = hgvs_p or ""
    if "fs" in hp or "frameshift" in hp.lower():
        return "frameshift_variant"
    if re.search(r"(Ter|\*)\)?$", hp) or hp.endswith("*"):
        return "stop_gained"
    if "Met1" in hp or "M1" in hp:
        return "start_lost"
    if re.search(r"del|dup|ins", hp):
        return "inframe_indel" if abs(len(ref) - len(alt)) % 3 == 0 and len(ref) != len(alt) else "frameshift_variant"
    return "missense_variant"


def _used(repo: Path):
    """Genes, coordinates and background samples already committed to the answer key."""
    genes, coords, samples = set(), set(), set()
    pv = repo / "manifest" / "planted_variants.tsv"
    if pv.exists():
        for r in csv.DictReader(open(pv), delimiter="\t"):
            genes.add((r["gene"] or "").upper())
            coords.add((_norm(r["chrom"]), str(r["pos"]), r["ref"].upper(), r["alt"].upper()))
    coh = repo / "manifest" / "cohort.tsv"
    if coh.exists():
        for r in csv.DictReader(open(coh), delimiter="\t"):
            genes.add((r["gene"] or "").upper())
            samples.add(r["sample"])
    return genes, coords, samples


def _next_syn_id(repo: Path) -> int:
    mx = 0
    coh = repo / "manifest" / "cohort.tsv"
    if coh.exists():
        for r in csv.DictReader(open(coh), delimiter="\t"):
            m = re.match(r"SYN-(\d+)", r["syn_id"] or "")
            if m:
                mx = max(mx, int(m.group(1)))
    return mx + 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenopackets", required=True, help="dir of GA4GH phenopacket *.json (recursive)")
    ap.add_argument("--backgrounds", required=True, help="file of candidate 1000G sample ids (one/line)")
    ap.add_argument("--n", type=int, default=100, help="how many new cases to select")
    ap.add_argument("--min-hpo", type=int, default=3, help="require at least this many HPO terms")
    ap.add_argument("--start", type=int, default=0, help="first new SYN id number (0 = max existing + 1)")
    ap.add_argument("--repo", default=str(REPO), help="repo root (to read the existing answer key)")
    ap.add_argument("--out-cohort", required=True)
    ap.add_argument("--out-plan", required=True)
    a = ap.parse_args(argv)

    repo = Path(a.repo)
    used_genes, used_coords, used_samples = _used(repo)
    backgrounds = [b for b in Path(a.backgrounds).read_text().split() if b and b not in used_samples]
    start = a.start or _next_syn_id(repo)

    files = sorted(glob.glob(os.path.join(a.phenopackets, "**", "*.json"), recursive=True))
    print(f"scanning {len(files)} phenopackets; {len(used_genes)} genes / {len(used_coords)} coords already used",
          file=sys.stderr)

    picked = []                                        # one dict per selected case
    seen_gene: set[str] = set()
    for fp in files:
        if len(picked) >= a.n:
            break
        try:
            d = load_phenopacket(fp)
        except Exception:
            continue
        hpo = d["hpo_terms"]
        if len(hpo) < a.min_hpo or not d["disease"]["label"]:
            continue
        caus = [v for v in d["variants"] if v["status"] == "CAUSATIVE" and v["gene"]]
        # primary = first causative variant on a fresh gene + coordinate
        prim = None
        for v in caus:
            g = v["gene"].upper()
            key = (_norm(v["chrom"]), str(v["pos"]), v["ref"], v["alt"])
            if g in used_genes or g in seen_gene or key in used_coords:
                continue
            prim = v
            break
        if not prim:
            continue
        seen_gene.add(prim["gene"].upper())
        others = [v for v in caus if (v["chrom"], v["pos"], v["ref"], v["alt"])
                  != (prim["chrom"], prim["pos"], prim["ref"], prim["alt"])]
        picked.append({"pkt": os.path.relpath(fp, a.phenopackets), "disease": d["disease"]["label"],
                       "hpo": hpo, "prim": prim, "others": others})

    if len(picked) < a.n:
        print(f"WARN: only {len(picked)} fresh cases found (asked for {a.n}); "
              f"widen --phenopackets or lower --min-hpo.", file=sys.stderr)
    if len(picked) > len(backgrounds):
        print(f"ERROR: {len(picked)} cases but only {len(backgrounds)} unused backgrounds — "
              f"add more sample ids to {a.backgrounds}.", file=sys.stderr)
        return 1

    cols = ["syn_id", "sample", "gene", "chrom", "pos", "ref", "alt", "consequence", "disease", "hpo"]
    plan = {}
    with open(a.out_cohort, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for i, c in enumerate(picked):
            sid = f"SYN-{start + i:03d}"
            p = c["prim"]
            w.writerow([sid, backgrounds[i], p["gene"], "chr" + _norm(p["chrom"]), p["pos"],
                        p["ref"], p["alt"], _coarse_consequence(p["hgvs_p"], p["ref"], p["alt"]),
                        c["disease"], ",".join(c["hpo"])])
            # faithful genotype
            if (p["zygosity"] or "").startswith("hom"):
                plan[sid] = {"mode": "hom"}
            elif c["others"]:
                o = c["others"][0]
                plan[sid] = {"mode": "compound_het", "chrom": "chr" + _norm(o["chrom"]),
                             "pos2": o["pos"], "ref2": o["ref"], "alt2": o["alt"],
                             "zyg2": o.get("zygosity") or "het"}
            else:
                plan[sid] = {"mode": "single_het"}
    json.dump(plan, open(a.out_plan, "w"), indent=0)

    from collections import Counter
    modes = Counter(v["mode"] for v in plan.values())
    print(f"selected {len(picked)} -> {a.out_cohort} (SYN-{start:03d}..SYN-{start + len(picked) - 1:03d})",
          file=sys.stderr)
    print(f"  genotypes: {dict(modes)}", file=sys.stderr)
    print(f"  NB: run annotate + reconcile_consequence.py to set the authoritative consequence.",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
