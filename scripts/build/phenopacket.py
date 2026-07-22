#!/usr/bin/env python3
"""Read a GA4GH Phenopacket (v2 JSON) into the fields the cohort builder needs.

Self-contained (standard library only) — a corrected fork of vcf2report's loader.

Extracts, per case:
  - subject_id
  - hpo_terms   : the case's FULL, ordered HPO set (excluded features skipped)
  - disease     : {"label", "id"} — the diagnosis (interpretations[].diagnosis.disease,
                  falling back to a top-level diseases[] term)
  - variants    : list of causative variants (GRCh38 vcfRecord coord + gene + HGVS +
                  zygosity from allelicState + interpretationStatus)

CORRECTIONS baked in (see scripts/build/BUILD.md → "Pitfalls"):
  * disease IS extracted (the original loader returned only hpo/variants, so every
    expansion case shipped an empty disease column).
  * the FULL HPO set is returned — never a 6-term subset.
  * disease/HPO are returned verbatim — never truncated.
  * only GRCh38 vcfRecords are emitted, and `interpretationStatus` is surfaced so the
    builder can require a CAUSATIVE variant.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_ZYGOSITY = {
    "GENO:0000135": "het", "heterozygous": "het",
    "GENO:0000136": "hom", "homozygous": "hom",
    "GENO:0000134": "hemi", "hemizygous": "hemi",
}
# GRCh38 assembly aliases; a vcfRecord naming another build is skipped so a GRCh37
# coordinate can never collide with a GRCh38 one.
_GRCH38 = {"grch38", "hg38", "grch38.p13", "grch38.p14", "grch38.p12", "grch38.p7", ""}


def _hpo_terms(pkt: dict) -> list[str]:
    terms: list[str] = []
    for feat in pkt.get("phenotypicFeatures", []) or []:
        if feat.get("excluded"):
            continue
        tid = (feat.get("type") or {}).get("id")
        if tid and tid.startswith("HP:"):
            terms.append(tid)
    return terms


def _disease(pkt: dict) -> dict:
    """The case diagnosis {label, id}. Prefer an interpretation's diagnosis.disease;
    fall back to a top-level diseases[] term. Empty label if none is recorded."""
    for interp in pkt.get("interpretations", []) or []:
        dz = (interp.get("diagnosis") or {}).get("disease") or {}
        if dz.get("label"):
            return {"label": dz["label"].strip(), "id": dz.get("id", "")}
    for dz in pkt.get("diseases", []) or []:
        if dz.get("excluded"):
            continue
        term = dz.get("term") or {}
        if term.get("label"):
            return {"label": term["label"].strip(), "id": term.get("id", "")}
    return {"label": "", "id": ""}


def _variant(gi: dict) -> dict | None:
    vi = gi.get("variantInterpretation") or {}
    desc = vi.get("variationDescriptor") or {}
    rec = desc.get("vcfRecord") or {}
    if str(rec.get("genomeAssembly", "")).lower() not in _GRCH38:
        return None
    if not rec.get("chrom") or rec.get("pos") in (None, "") or not rec.get("ref") or not rec.get("alt"):
        return None
    exprs = {e.get("syntax"): e.get("value") for e in desc.get("expressions", []) or []}
    state = desc.get("allelicState") or {}
    return {
        "chrom": str(rec["chrom"]),
        "pos": int(rec["pos"]),
        "ref": str(rec["ref"]).upper(),
        "alt": str(rec["alt"]).upper(),
        "gene": (desc.get("geneContext") or {}).get("symbol"),
        "hgvs_c": exprs.get("hgvs.c"),
        "hgvs_p": exprs.get("hgvs.p"),
        "zygosity": _ZYGOSITY.get(state.get("id")) or _ZYGOSITY.get((state.get("label") or "").lower()),
        "status": gi.get("interpretationStatus", ""),   # CAUSATIVE / CANDIDATE / ...
    }


def load_phenopacket(path: str | Path) -> dict[str, Any]:
    """Return {'subject_id', 'hpo_terms', 'disease', 'variants'}."""
    pkt = json.loads(Path(path).read_text())
    variants = []
    for interp in pkt.get("interpretations", []) or []:
        for gi in (interp.get("diagnosis") or {}).get("genomicInterpretations", []) or []:
            v = _variant(gi)
            if v:
                variants.append(v)
    return {
        "subject_id": (pkt.get("subject") or {}).get("id") or "PHENOPACKET",
        "hpo_terms": _hpo_terms(pkt),
        "disease": _disease(pkt),
        "variants": variants,
    }


if __name__ == "__main__":
    import sys
    d = load_phenopacket(sys.argv[1])
    print(f"subject : {d['subject_id']}")
    print(f"disease : {d['disease']['label']}  ({d['disease']['id']})")
    print(f"hpo ({len(d['hpo_terms'])}): {','.join(d['hpo_terms'])}")
    for v in d["variants"]:
        print(f"  {v['gene']}  {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
              f"{v['zygosity']}  {v['status']}")
