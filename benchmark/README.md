# Benchmark — reference results

**These files are for validating variant-interpretation tools.** They record the **expected result of
one specific tool** — the [vcf2report](https://github.com/gbbarra/vcf2report) ACMG engine (a *separate*
repository) — run over this cohort. vcf2report uses them as its own validation set; anyone scoring a
different annotator / classifier can compare against them, but should expect different numbers.

> ⚠️ **Not part of the answer key.** The dataset's ground truth — gene · coordinate · HPO · disease,
> in `manifest/` + `sidecars/` — is **tool-agnostic**. The per-case **`outcome` / `tier`** here is
> whatever *vcf2report* produced; a different ACMG engine will classify the same variant differently.
> Score your tool against the **manifest** (did it recover the planted gene?), not against these tiers.

## Files

| File | What it is |
|---|---|
| `vcf2report.data-v1.tsv` | per-case result of the vcf2report engine on the `data-v1` cohort (200 rows) |

## Result (vcf2report on `data-v1`)

Measured on the current cohort (full-HPO sidecars + SnpEff-reconciled consequence): the planted variant
reaches the diagnostic (**primary**) finding in **177 / 200 (88.5%)**.

| `outcome` | n | meaning |
|---|---:|---|
| `primary` | **177** | the planted gene is the top diagnostic finding |
| `other` | 10 | planted gene reported, but below the top finding |
| `probable_vus` | 5 | held at VUS, flagged probably-pathogenic |
| `absent` | 7 | planted gene not surfaced |
| `carrier` | 1 | one allele of a recessive case (SPINT2) |

Of the 177 `primary`, vcf2report's ACMG tier is **85 Pathogenic · 22 Likely Pathogenic · 70 VUS** —
i.e. the gene is recovered as the lead finding even when the engine conservatively holds the tier at VUS
(the cohort deliberately oversamples missense / in-frame, which produce VUS without corroboration).

### Why the 23 non-`primary` cases are not misses

Each of the 23 was cross-checked against the **manifest** (the real SnpEff consequence + ClinVar
significance). Every one is a known, disclosed limitation — **none is a regression**; they split
cleanly by *why* the engine held back:

| Category | n | Cases | consequence · ClinVar | Why it is expected |
|---|---:|---|---|---|
| Non-coding RNA | 2 | RNU5B-1, RNU4-2 | `non_coding_transcript_exon_variant` | outside ACMG's protein-coding scope |
| Low-impact splice-region | 3 | ADA, BBS1, BNIP1 | `splice_region_variant` | the real consequence is splice-**region** (not donor/acceptor); the impact filter doesn't promote it |
| ClinVar-benign plant | 1 | ZIC2 | `Benign/Likely_benign` | ClinVar calls it benign — the engine is *right* not to call it |
| Under-characterized gene | 1 | C10ORF71 | frameshift · Pathogenic | an ORF with no HPO/constraint coverage, so PVS1 can't engage — a **data gap**, not an error |
| Carrier (recessive) | 1 | SPINT2 | frameshift · one allele | a lone het in a recessive gene is a *carrier*, not a diagnosis |
| Missense / in-frame at VUS | 7 | U2AF2, SETD2, TRAF7, MAPK8IP3, MOCS3, ETFB, AMT | missense / in-frame | held at VUS without corroborating evidence |
| Likely-Pathogenic, off-phenotype | 2 | ANTXR2, BRPF1 | frameshift · Likely Pathogenic | reached LP but HPO overlap fell below the primary threshold → surfaced as `other` |
| Probable-VUS (triaged) | 5 | RBSN, ELFN1, ATP6V1A, SPG7, CCIN | missense / in-frame | VUS carrying a molecular signal, flagged for review |

So **177/200** are recovered as the lead finding, and the remaining 23 are non-coding / low-impact-splice
/ ClinVar-benign / data-gap / carrier / held-at-VUS — all disclosed, none a defect in the cohort. (This
breakdown was produced with the vcf2report engine; see its `docs/BENCHMARK.md` for the full harness.)

## Column dictionary — `vcf2report.data-v1.tsv`

| Column | Meaning |
|---|---|
| `syn_id` · `gene` | the case and its planted gene |
| `outcome` | where the planted gene landed: `primary` \| `other` \| `probable_vus` \| `carrier` \| `absent` |
| `tier` | vcf2report's ACMG classification of the planted variant (blank when `absent`) |
| `candidates` | number of candidate variants the engine considered in that exome |
| `error` | run error, if any (blank = clean) |

## Reproduce

Run vcf2report on each `realistic/SYN-NNN.vcf.gz` + `sidecars/SYN-NNN.hpo.txt`, then compare the top
finding's gene to `manifest/planted_variants.tsv`. See the [vcf2report](https://github.com/gbbarra/vcf2report)
repo for the engine and its harness. Re-measure after any cohort change — these numbers are
cohort- and tool-specific.
