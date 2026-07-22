# hpo-spiked-exomes

**A ground-truth benchmark for variant annotators, variant callers, and clinical-report generators.**

200 **real, healthy [1000 Genomes](https://www.internationalgenome.org/) exomes**, each with **one
known pathogenic variant** (plus its real second allele where the source case is biallelic) **planted
at an exact coordinate**, carrying that case's **HPO phenotype**. The answer key is known — but in the
**realistic** build the plant is **not marked** in the VCF, so any tool can be scored **blind**.

> ⚠️ **Synthetic, not real patients.** Real public backgrounds with an inserted variant. De-identified.
> **Not for clinical use.**

---

## Get the data

The VCFs (~3 GB) ship as **release assets**; this git repo holds the docs, the answer key, and the
scripts.

```bash
git clone https://github.com/gbbarra/hpo-spiked-exomes.git && cd hpo-spiked-exomes
bash fetch.sh              # downloads + checksums + extracts the release assets into ./realistic and ./realistic_annotated
```

After `fetch.sh`:

| Path | What it is |
|---|---|
| `realistic/SYN-NNN.vcf.gz` | **Raw, tell-free** VCF — the plant carries a real DRAGEN call's INFO/FORMAT and **no marker**. Un-annotated. |
| `realistic_annotated/SYN-NNN.annotated.vcf.gz` | The **same** VCF **SnpEff-annotated** (`GRCh38.mane.1.5.refseq`) — `ANN/LOF/NMD` added, still tell-free. |
| `manifest/planted_variants.tsv` | **The answer key** — every planted allele: `chrom:pos:ref:alt`, gene, zygosity, primary/second, consequence, ClinVar significance/review/id. |
| `manifest/cohort.tsv` | Per-sample config (sample id, gene, coord, consequence, disease, HPO). |
| `sidecars/SYN-NNN.planted.tsv` · `sidecars/SYN-NNN.hpo.txt` | Per-sample answer key + HPO terms (one `HP:` per line). |

The raw VCFs are **byte-derivable** from the annotated (`bcftools annotate -x INFO/ANN,INFO/LOF,INFO/NMD`),
so the two forms are guaranteed consistent.

## What's inside: 200 cases

- **Backgrounds** — 200 **distinct** 1000 Genomes **DRAGEN v4.4.7** exomes (public AWS Open Data bucket
  `1000genomes-dragen-v4-4-7`), streamed, normalized, and **subset to the MANE/GENCODE exome BED**
  (~100k variants each). Every case uses a **different** sample across diverse populations; no
  background is reused.
- **Plants** — each causative variant is the **exact `chrom:pos:ref:alt`** from a real
  **[GA4GH Phenopacket Store](https://github.com/monarch-initiative/phenopacket-store) 0.1.27** case,
  with that case's **HPO terms**, gene, consequence, and disease. Real ClinVar `CLNSIG` where the
  coordinate is in ClinVar; a **synthetic** label (flagged in the manifest) where it is not.
- **Faithful genotypes** — the patient's real zygosity: **40 compound-heterozygous** (both true
  alleles), **75 homozygous**, **85 single-allele**.
- **Consequence spread** — stratified across missense / stop-gained / frameshift / in-frame / start-loss.

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
- Some plants carry a **synthetic ClinVar label** (coordinate absent from the ClinVar release) — flagged
  in `manifest/planted_variants.tsv` (`clnrevstat` / no real assertion).
- **GRCh38 only.**

## Sources & citations

- **1000 Genomes / IGSR**, Illumina DRAGEN v4.4.7 re-analysis (`1000genomes-dragen-v4-4-7`, AWS Open Data).
- **GA4GH Phenopacket Store** (Danis et al., 2023) — the causative variants + HPO.
- **ClinVar** (NCBI, public domain) · **gnomAD v4.1** · **MANE / GENCODE** · **Human Phenotype Ontology**.
- Annotation: **SnpEff 5.4c** (`GRCh38.mane.1.5.refseq`).

Reproducibility scripts under `scripts/` (the tell-free transform + HPO sidecar generator); the full
build pipeline is in [vcf2report](https://github.com/gbbarra/vcf2report).

## License

Tooling and the answer key: **MIT** (see `LICENSE`). The VCF data derives from public sources, each
under its own license (1000 Genomes open; Phenopacket Store open; ClinVar public domain). No
AlphaMissense data is redistributed here.
