# Building more cases — the cohort expansion pipeline

Tooling to grow the benchmark: take real **1000G DRAGEN** exomes + real **GA4GH Phenopacket-Store**
cases and produce more tell-free planted VCFs **plus** a correct answer key. Imported from the
`vcf2report` build tree and corrected so the mistakes we hit building `data-v1` cannot recur — every
one is listed under [Pitfalls](#pitfalls--do-not-regress) with the guard that prevents it.

The annotation mechanism is **SnpEff** (`GRCh38.mane.1.5.refseq`) throughout. The answer key's
`consequence` is *defined* as what SnpEff calls at the planted coordinate — never a hand-written or
heuristic value (that was a real bug; see pitfall #4).

## Prerequisites (not in this repo)

| Need | Where |
|---|---|
| Phenopacket Store 0.1.27 | `git clone --depth 1 https://github.com/monarch-initiative/phenopacket-store` |
| MANE/GENCODE exome BED | `exome_hg38.bed` (vcf2report `data/gnomad/`) — the region the engine covers |
| ClinVar GRCh38 VCF | `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` |
| SnpEff + `GRCh38.mane.1.5.refseq` | `$SNPEFF_JAR` or `snpEff` on PATH — **same DB as the shipped data** |
| CLI | `curl bcftools bgzip tabix zstd python3 java gh` |
| backgrounds.txt | candidate 1000G sample ids (one per line); unused ones are assigned |

## The pipeline (add N cases, e.g. SYN-201…)

```bash
# 1. pick fresh cases — stores disease + FULL HPO, requires a CAUSATIVE variant,
#    excludes every gene/coord/background already in the answer key.
python3 scripts/build/select_cases.py \
    --phenopackets ~/phenopacket-store --backgrounds backgrounds.txt --n 100 \
    --out-cohort new_cohort.tsv --out-plan new_plan.json

# 2. plant the exact variant into each 1000G exome (marker-bearing, QC-passing template)
BED=/path/exome_hg38.bed CLINVAR_VCF=/path/clinvar.vcf.gz \
    bash scripts/build/make_cohort.sh new_cohort.tsv raw/

# 3. apply the FAITHFUL genotype (compound-het / hom / single-het) from the plan
python3 scripts/build/build_biallelic.py \
    --plan new_plan.json --cohort-tsv new_cohort.tsv --src-dir raw/ --out biallelic/

# 4. make each record tell-free + capture the truth into the manifest & per-case sidecar
for v in biallelic/SYN-*.synthetic.vcf.gz; do sid=$(basename "$v" .synthetic.vcf.gz)
  python3 scripts/realisticize_cohort.py "$v" realistic/$sid.vcf \
      --syn-id "$sid" --sidecar sidecars/$sid.planted.tsv --manifest manifest/planted_variants.tsv
  bgzip -f realistic/$sid.vcf; done
tail -n +2 new_cohort.tsv >> manifest/cohort.tsv          # append the new cohort rows

# 5. annotate with SnpEff (the benchmark's annotator)
for v in realistic/SYN-*.vcf.gz; do sid=$(basename "$v" .vcf.gz)
  bash scripts/build/annotate.sh "$v" realistic_annotated/$sid.annotated.vcf.gz; done

# 6. set the AUTHORITATIVE consequence from the real annotation (cohort + manifest + sidecars)
python3 scripts/reconcile_consequence.py --annotated realistic_annotated

# 7. materialise the HPO sidecars (one HP: per line, matching cohort.tsv)
python3 scripts/fill_hpo_sidecars.py manifest/cohort.tsv sidecars

# 8. pack clean tarballs + checksums (+ upload)
bash scripts/build/publish.sh realistic/ realistic_annotated/ data-v2 --upload
```

Then re-run the answer-key audit (the same checks used in review) and the benchmark before announcing.

## Files

| Stage | Script | Notes |
|---|---|---|
| load | `build/phenopacket.py` | self-contained loader — disease + full HPO + causative variants |
| select | `build/select_cases.py` | fresh-case picker, no reuse, faithful-genotype plan |
| plant | `../spike_variant.py` (+`spike_pathogenic.py`) | exact coord, QC-passing borrowed template |
| genotype | `build/build_biallelic.py` | compound-het / hom from the plan |
| tell-free | `../realisticize_cohort.py` | strips markers, writes truth sidecar + manifest |
| annotate | `build/annotate.sh` | bcftools norm + SnpEff (chr-naming guard) |
| consequence | `../reconcile_consequence.py` | sets the answer key to the real SnpEff term |
| sidecars | `../fill_hpo_sidecars.py` | HPO sidecars from cohort.tsv |
| publish | `build/publish.sh` | clean tarballs (`COPYFILE_DISABLE`) + `shasum -a 256` |

## Pitfalls — do NOT regress

Every item below was a real defect found while building/reviewing `data-v1`. The guard is in place; keep it.

1. **Empty `disease`.** The old loader never read the diagnosis and the selector had dead code, so 100
   cases shipped blank. → `phenopacket._disease()` reads `interpretations[].diagnosis.disease`;
   `select_cases` stores it. Never ship a blank disease.
2. **HPO capped at 6 terms.** → `select_cases` stores the **full, ordered** HPO set; the `.hpo.txt`
   sidecar must equal the `cohort.tsv` `hpo` column (audit checks this).
3. **Disease truncated at 60 chars.** → written verbatim everywhere; never `[:60]` a label.
4. **`consequence` = a HGVS heuristic.** The selector's guess disagreed with the real annotator in
   ~34 cases. → the value in `cohort.tsv`/manifest is **only a placeholder** until
   `reconcile_consequence.py` overwrites it from SnpEff. cohort.tsv and the manifest must AGREE
   (reconcile touches both). Score against the SnpEff term, never the placeholder.
5. **`clnrevstat` encoding.** Synthetic labels once used spaces (`criteria provided, single submitter`)
   vs ClinVar's underscores. → `spike_variant.py` emits the underscored form.
6. **Non-causative variant.** → `select_cases` requires `interpretationStatus == CAUSATIVE`.
7. **Reusing a gene / coordinate / background.** → `select_cases` excludes everything already in
   `manifest/planted_variants.tsv` + `cohort.tsv`; new SYN ids continue after the max.
8. **QC-failing borrowed template.** A low-GQ template silently dropped 17 cases at the engine's QC.
   → `spike_variant._pick_template` requires **GQ ≥ 30, DP ≥ 25, balanced het**.
9. **Markers leaking into the VCF.** → `realisticize_cohort.py` strips all `SPIKED/GENE/CSQ/CLN*`;
   truth lives ONLY in `manifest/` + `sidecars/`, never in the shipped VCF.
10. **Hand-written consequence on the 2nd allele.** → `build_biallelic` writes only `GENE`+`SPIKED2`
    for the trans allele, so it can't disagree with SnpEff; its consequence is re-derived on annotate.
11. **Empty / mismatched HPO sidecar.** An empty `--hpo` file yields a silent genotype-only run
    (no PP4). → `fill_hpo_sidecars.py`; sidecar must equal the cohort HPO.
12. **macOS `._*` AppleDouble junk in tarballs.** `--exclude='._*'` does NOT stop bsdtar synthesizing
    them from xattrs. → `publish.sh` sets `COPYFILE_DISABLE=1` **and** verifies 0 `._*` before upload.
13. **Silent annotation failure (chromosome naming).** A UCSC/Ensembl mismatch makes SnpEff annotate
    nothing (every record `ERROR_CHROMOSOME_NOT_FOUND`). → `annotate.sh` renames in/out and **fails
    loudly** if fewer than half the records get `ANN`.
14. **Checksum format.** `fetch.sh` verifies with `shasum -a 256 -c` / `sha256sum -c`. → `publish.sh`
    writes `SHA256SUMS` with `shasum -a 256`.
15. **SnpEff DB drift.** A different DB across releases makes consequences incomparable. → keep
    `SNPEFF_DB=GRCh38.mane.1.5.refseq` (the shipped DB) for every expansion.

## Versioning

Publish an expansion either as a **new tag** (`data-v2`, additive) or by re-packing the combined
cohort into `data-v1`. Either way: bump the README's cohort counts, the ClinVar-label distribution,
the composition table, and re-measure the reference recovery — those numbers are cohort-specific.
