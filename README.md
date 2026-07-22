# hpo-spiked-exomes

[![Release](https://img.shields.io/github/v/release/gbbarra/hpo-spiked-exomes?label=data&color=blue)](https://github.com/gbbarra/hpo-spiked-exomes/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
![Assembly: GRCh38](https://img.shields.io/badge/assembly-GRCh38-informational)
![Cases: 200](https://img.shields.io/badge/cases-200-blueviolet)

**A ground-truth benchmark for variant annotators, variant callers, and clinical-report generators.**

200 **real, healthy [1000 Genomes](https://www.internationalgenome.org/) exomes**, each carrying **the
causative variant of a real, published, HPO-phenotyped case** (plus its real second allele where the
source case is biallelic) **planted at an exact coordinate**. The answer key is known — but in the
**realistic** build the plant is **not marked** in the VCF, so any tool can be scored **blind**.

> ℹ️ **What "causative" means here.** The plant is the exact variant the source Phenopacket case was
> diagnosed on — not necessarily a variant ClinVar currently calls Pathogenic. Most plants are ClinVar
> P/LP, but **~29 carry a real ClinVar label that is VUS / Conflicting / Benign** and **32 have no
> ClinVar record at all** (synthetic label). The ground truth is *gene + coordinate + phenotype*, not
> a ClinVar lookup — see [The answer key](#the-answer-key). This is deliberate: it tests
> phenotype-driven interpretation, not ClinVar recall.

> ⚠️ **Synthetic, not real patients.** Real public backgrounds with an inserted variant. De-identified.
> **Not for clinical use.**

---

## Get the data

The VCFs (~3 GB) ship as **[release assets](https://github.com/gbbarra/hpo-spiked-exomes/releases)**;
this git repo holds the docs, the answer key, and the scripts.

```bash
git clone https://github.com/gbbarra/hpo-spiked-exomes.git && cd hpo-spiked-exomes
bash fetch.sh              # downloads + checksums + extracts the release assets
# bash fetch.sh data-v1   # a specific release tag (default: data-v1)
```

`fetch.sh` needs **`zstd`** and a checksum tool (**`shasum`** or **`sha256sum`**), and uses the GitHub
CLI (`gh`) when present, else `curl`:

| OS | Install prerequisites |
|---|---|
| Debian / Ubuntu | `apt-get install zstd`  (`sha256sum` ships with coreutils) |
| macOS (Homebrew) | `brew install zstd`  (`shasum` ships with the system perl) |

After `fetch.sh`:

| Path | What it is |
|---|---|
| `realistic/SYN-NNN.vcf.gz` | **Raw, tell-free** VCF — the plant carries a real DRAGEN call's INFO/FORMAT and **no marker**. Un-annotated. |
| `realistic_annotated/SYN-NNN.annotated.vcf.gz` | The **same** VCF **SnpEff-annotated** (`GRCh38.mane.1.5.refseq`) — `ANN/LOF/NMD` added, still tell-free. |
| `manifest/planted_variants.tsv` | **The answer key** — every planted allele (see the schema below). |
| `manifest/cohort.tsv` | Per-sample config (sample id, gene, coord, consequence, disease, HPO). |
| `sidecars/SYN-NNN.planted.tsv` · `sidecars/SYN-NNN.hpo.txt` | Per-sample answer key + HPO terms (one `HP:` per line). |

The raw VCFs are **byte-derivable** from the annotated (`bcftools annotate -x INFO/ANN,INFO/LOF,INFO/NMD`),
so the two forms are guaranteed consistent.

## Quick start — validate on one case

```bash
# SYN-004 — NIPBL / Cornelia de Lange syndrome 1
your-report-pipeline \
    --vcf realistic/SYN-004.vcf.gz \
    --hpo sidecars/SYN-004.hpo.txt
# expected top finding:  NIPBL  chr5:37022325 C>T  stop_gained
#   (manifest: Pathogenic, "Cornelia de Lange syndrome 1")
```

Then score across all 200 with `manifest/planted_variants.tsv` as the key. **Score blind on
`realistic/`** — the plant carries no synthetic hint.

## What's inside: 200 cases

- **Backgrounds** — 200 **distinct** 1000 Genomes **DRAGEN v4.4.7** exomes (public AWS Open Data bucket
  `1000genomes-dragen-v4-4-7`), streamed, normalized, and **subset to the MANE/GENCODE exome BED**
  (~100k variants each). Every case uses a **different** sample across diverse populations; no
  background is reused.
- **Plants** — each causative variant is the **exact `chrom:pos:ref:alt`** from a real
  **[GA4GH Phenopacket Store](https://github.com/monarch-initiative/phenopacket-store) 0.1.27** case,
  with that case's **HPO terms**, gene, consequence, and disease. Real ClinVar `CLNSIG` where the
  coordinate is in ClinVar; a **synthetic** label (flagged in the manifest) where it is not.
- **Faithful genotypes** — the patient's real zygosity.
- **Consequence spread** — stratified across missense / stop-gained / frameshift / in-frame / start-loss.

### Cohort composition

| Molecular consequence | n | | Genotype of the plant | n |
|---|---:|---|---|---:|
| missense | 91 | | single-allele heterozygous | 85 |
| frameshift | 41 | | **compound heterozygous** (both true alleles) | 40 |
| stop-gained | 38 | | homozygous | 75 |
| in-frame indel | 21 | | | |
| start-loss | 9 | | | |
| **Total primary plants** | **200** | | **Total** | **200** |

240 planted alleles in all = 200 primary + 40 trans (`second`) alleles for the compound-het cases.

## The answer key

`manifest/planted_variants.tsv` is the source of truth; `sidecars/SYN-NNN.planted.tsv` is the per-case slice.

| Column | Meaning | Values |
|---|---|---|
| `syn_id` | Case id | `SYN-001` … `SYN-200` |
| `chrom` `pos` `ref` `alt` | GRCh38 coordinate of the planted allele | e.g. `chr5` `37022325` `C` `T` |
| `gene` | Gene symbol | HGNC symbol |
| `zygosity` | Genotype of the planted call | `het` \| `hom` |
| `allele` | Role of this row | `primary` (the reported causative allele) \| `second` (trans allele of a compound-het case) |
| `consequence` | Molecular consequence *(primary rows only)* | `missense_variant`, `stop_gained`, `frameshift_variant`, `inframe_indel`, `start_lost`, … |
| `clnsig` | ClinVar significance carried from the exact coordinate; **synthetic `Pathogenic`** when `clnvid` is empty | ClinVar CLNSIG string |
| `clnrevstat` | ClinVar review status *(underscored ClinVar encoding)* | e.g. `criteria_provided,_multiple_submitters,_no_conflicts` |
| `clnvid` | ClinVar Variation ID; **empty ⇒ coordinate absent from the ClinVar release ⇒ `clnsig` is synthetic** | integer or empty |

`second`-allele rows leave `consequence` / `clnsig` / `clnrevstat` / `clnvid` blank.

### ClinVar label of the primary plant

Because the plant is the *phenopacket* causative variant (not a ClinVar pick), the carried ClinVar label
varies. Filter on `clnvid` to separate real assertions from synthetic ones:

| `clnsig` of the primary plant | n | Note |
|---|---:|---|
| Pathogenic | 135 | **32 of these are synthetic** (coordinate absent from ClinVar → default `Pathogenic`, `clnvid` empty) |
| Pathogenic/Likely_pathogenic | 18 | |
| Likely_pathogenic | 18 | |
| Uncertain_significance | 13 | real ClinVar VUS |
| Conflicting_classifications_of_pathogenicity | 12 | real ClinVar conflict |
| Benign/Likely_benign · risk_factor · drug_response · not_provided | 4 | real ClinVar, non-pathogenic label |

So **139 plants carry a real ClinVar P/LP assertion**, **32 carry a synthetic label**, and **29 carry a
real ClinVar label that is not P/LP** — expected, since the truth is the phenotype-linked causative
variant, not ClinVar's current call.

## The realistic (tell-free) transform — and the honest adjustments

The naive way to spike a variant leaves obvious markers (`SPIKED=1`, `GENE/CSQ/CLN*` INFO, a minimal
`GT:DP:GQ:AD` FORMAT). The **realistic** build removes all of that: it **borrows a real background
call's full DRAGEN INFO/FORMAT** of the same zygosity, relocates it to the plant's coordinate/alleles,
and strips every marker. **The answer key lives only in this repo** (`manifest/`, `sidecars/`), never
in the VCF.

Two adjustments were needed, disclosed transparently:

1. **QC-passing templates.** A borrowed real call can have a low genotype quality; the plant would then
   inherit it and be dropped at a caller/engine's QC. The template is required to pass QC comfortably
   (**GQ ≥ 30, DP ≥ 25, balanced het**). *(Found when a first build silently lost 17 cases.)*
2. **De-circularization — a feature.** Because the tell-free VCF carries **no** synthetic `CSQ`/`CLNSIG`,
   a tool must work from the **real sequence context**, not from the spike's hints. A few plants whose
   marker-bearing classification leaned on a synthetic tag (a `CSQ=missense` the real annotator scores
   as a low-impact splice-region variant; a `CLNSIG` for a variant absent from real ClinVar) land on
   their honest tier here — a **cleaner, less circular** benchmark.

## How to validate your tool

- **Annotator** — annotate `realistic/*.vcf.gz`; for each planted coordinate (manifest) check the gene /
  consequence / HGVS you produce.
- **Caller** — the plant is a genuine-looking call; check it is recovered.
- **Report generator** — run your pipeline on a sample + its `sidecars/SYN-NNN.hpo.txt`; the expected
  finding is the planted gene (manifest). Score blind on `realistic/`.

Reference: with the [vcf2report](https://github.com/gbbarra/vcf2report) ACMG engine, the planted variant
reaches the diagnostic (primary) finding in **178/200** cases; the rest are honest limitations
(non-coding-RNA plants, HPO-unlinked genes, sub-threshold phenotype, missense held at VUS without
corroboration) — documented, not hidden.

## Honest limitations

- Synthetic: real backgrounds, **inserted** causative variant — tests classification/annotation, not
  real diagnostic yield.
- The plant is the **phenopacket** causative variant, so its ClinVar label is not uniform: **32/200**
  carry a **synthetic** ClinVar label (coordinate absent from the ClinVar release — `clnvid` empty) and
  **29/200** carry a **real but non-P/LP** label (VUS / Conflicting / Benign / …). Both are flagged in
  `manifest/planted_variants.tsv`.
- **GRCh38 only.**

## Reproducibility — versions

| Component | Version / source |
|---|---|
| Backgrounds | 1000 Genomes / IGSR, Illumina **DRAGEN v4.4.7** (`1000genomes-dragen-v4-4-7`, AWS Open Data) |
| Causative variants + HPO | **GA4GH Phenopacket Store 0.1.27** (Danis et al., 2023) |
| Functional annotation | **SnpEff 5.4c** (`GRCh38.mane.1.5.refseq`) |
| Allele frequencies | **gnomAD v4.1** |
| Transcripts / exome BED | **MANE / GENCODE** (1.5) |
| Clinical significance | **ClinVar** (NCBI, public domain, GRCh38) <!-- TODO: pin the exact dated ClinVar release used at build --> |
| Phenotype ontology | **Human Phenotype Ontology** |
| Assembly | **GRCh38** |

Reproducibility scripts under `scripts/`:

- `spike_variant.py` — plant one exact `chrom:pos:ref:alt` (the tell-free transform lives here).
- `spike_pathogenic.py` — plant *a* pathogenic ClinVar record for a gene (the marker-bearing base build).
- `realisticize_cohort.py` — strip the markers from an already-spiked VCF and emit the truth sidecar.
- `fill_hpo_sidecars.py` — materialise `SYN-NNN.hpo.txt` from the cohort TSV.

The full build pipeline is in [vcf2report](https://github.com/gbbarra/vcf2report).

## Citation

Barra, G. B. *hpo-spiked-exomes: a tell-free benchmark of HPO-linked pathogenic variants spiked into
real exomes* (2026). https://github.com/gbbarra/hpo-spiked-exomes

DOI: _pending Zenodo deposit._
<!-- Once minted, add the badge at the top: [![DOI](https://zenodo.org/badge/DOI/<doi>.svg)](https://doi.org/<doi>) -->

## License

Tooling and the answer key: **MIT** (see `LICENSE`) © Gustavo Barcelos Barra (gbbarra). The VCF data derives from public sources, each
under its own license (1000 Genomes open; Phenopacket Store open; ClinVar public domain). No
AlphaMissense data is redistributed here.
