<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-scalpel_logo_dark.png">
    <img alt="nf-core/scalpel" src="docs/images/nf-core-scalpel_logo_light.png">
  </picture>
</h1>

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/plasslab/nf-core-scalpel)

[![Get help on Slack](https://img.shields.io/badge/slack-nf--core%20%23scalpel-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/scalpel)
[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](https://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/scalpel** quantifies transcript isoforms and characterizes alternative polyadenylation (APA) at the single-cell level from 3'-tagged scRNA-seq data (10x Genomics Chromium 3', Drop-seq). It is the Nextflow (DSL2, nf-core template) implementation of [SCALPEL](https://github.com/plasslab/SCALPEL) ([Ake et al., _Nat Commun_ 2025](https://doi.org/10.1038/s41467-025-61118-0)).

From the FASTQ files and the Cell Ranger / STARsolo / Drop-seq outputs of each sample, the pipeline:

1. **Annotation processing** — quantifies transcripts in pseudobulk with salmon, averages the TPM over samples, and selects the isoforms distinguishable by 3'-tagged reads: isoforms with close 3' ends and identical exons near the 3' end are collapsed into the most expressed one.
2. **Reads processing** — splits the BAM files by chromosome (filtered cells, uniquely mapped reads, PCR duplicates removed), maps the reads to the selected isoforms (read agreement and distance filters) and removes internal priming artifacts. The internal priming reference is built from the genome by default.
3. **Isoform quantification** — learns the read distance distribution to the 3' end, computes fragment probabilities, estimates the isoform abundances of each cell and gene with an EM algorithm, and writes isoform-by-cell count matrices (iDGE).
4. **APA characterization** — differential isoform usage between samples or cell clusters, Seurat object and BAM files of the reads used.

<div align="center">
  <img src="docs/images/SCALPEL_diagram.png" alt="SCALPEL pipeline overview" width="800">
</div>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow (version >= 26.04). Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Quick test

```bash
# 10x Chromium test dataset (GSE104556, chr17:14-16Mb, a few minutes)
nextflow run plasslab/nf-core-scalpel -profile test,docker --outdir results_test

# same dataset in Drop-seq format
nextflow run plasslab/nf-core-scalpel -profile test_dropseq,docker --outdir results_test_dropseq
```

Use `conda`, `singularity` or `apptainer` instead of `docker` as needed. On Apple Silicon Macs with Docker, add `emulate_amd64` (e.g. `-profile test,docker,emulate_amd64`).

### Samplesheet

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cranger_path
SRR6129050,fastq/SRR6129050_R1.fastq.gz,fastq/SRR6129050_R2.fastq.gz,cellranger/SRR6129050
SRR6129051,fastq/SRR6129051_R1.fastq.gz,fastq/SRR6129051_R2.fastq.gz,cellranger/SRR6129051
```

| Column         | Description                                                                                                         |
| -------------- | ------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Sample name (no spaces).                                                                                            |
| `fastq_1`      | R1 FASTQ (`.fastq.gz`). Optional with `--isoform_annotation`.                                                       |
| `fastq_2`      | R2 FASTQ (`.fastq.gz`). Optional with `--isoform_annotation`.                                                       |
| `cranger_path` | Sample output directory: Cell Ranger, Cell Ranger ARC (multiome GEX) or STARsolo; Drop-seq: the BAM file.           |
| `dge`          | Drop-seq only: Drop-seq tools DGE count matrix (`DigitalExpression`).                                               |

- Paths can be absolute or relative to the samplesheet. The BAM, its index (created if missing) and the filtered cell barcodes are found automatically in the supported layouts.
- A sample sequenced in several runs or lanes is given as several rows with the same `sample` name and `cranger_path`: all its FASTQ pairs are quantified together.

### Running the pipeline

```bash
nextflow run plasslab/nf-core-scalpel \
   -profile <docker/singularity/conda/.../institute> \
   --samplesheet samplesheet.csv \
   --transcriptome gencode.vM21.transcripts.fa.gz \
   --gtf gencode.vM21.annotation.gtf.gz \
   --genome GRCm38.primary_assembly.genome.fa.gz \
   --outdir results
```

| Input             | Description                                                                                                                       |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `--transcriptome` | Transcript FASTA (e.g. GENCODE `transcripts.fa.gz`), **not** the genome.                                                          |
| `--gtf`           | Matching GTF annotation (same release).                                                                                           |
| `--genome`        | Genome FASTA, used to build the internal priming reference (A-rich genomic sites, saved to `<outdir>/ip_reference/` for reuse).   |
| `--ip_reference`  | Optional pre-built internal priming reference used instead of `--genome`: [mm10](https://zenodo.org/records/15664563/files/mm10_polya.track.tar.gz?download=1), [GRCh38](https://zenodo.org/records/15717592/files/hg38_ipriming_sites.bed.tar.gz?download=1), or `ip_reference.track.gz` of a previous run. |

> [!IMPORTANT]
> The BAM files, the GTF and the genome / internal priming reference must use the **same chromosome names** (e.g. `chr1` vs `1`). The pipeline checks this before processing the reads and stops with an explicit message otherwise.

Common use cases:

```bash
# compare cell types / clusters instead of samples (sample, barcode, cluster)
--clusters clusters.tsv

# restrict the analysis to selected cells (CSV: sample,barcodes)
--barcodes_whitelist barcodes_whitelist.csv

# add samples later on the same isoforms (no salmon, FASTQ optional: BAM input only)
--isoform_annotation results_run1/isoform --ip_reference results_run1/ip_reference/ip_reference.track.gz

# Drop-seq data
--sequencing_platform dropseq
```

### Parameters

| Parameter               | Default    | Original SCALPEL | Description                                                                                  |
| ----------------------- | ---------- | ---------------- | -------------------------------------------------------------------------------------------- |
| `--sequencing_platform` | `chromium` |                  | `chromium` (CB/UB tags) or `dropseq` (XC/XM tags).                                            |
| `--distance_3end`       | `100`      | `50`             | Maximum distance (nt) between the 3' ends of isoforms collapsed together.                    |
| `--distance_profile`    | `500`      | `600`            | Window (nt from the 3' end) in which collapsed isoforms must have identical exons.            |
| `--distance_reads`      | `600`      | `600`            | Maximum distance (nt) of a read from the isoform 3' end.                                      |
| `--min_mapq`            | `255`      | `0`              | Minimum mapping quality (255: uniquely mapped reads, as counted in the gene matrix).          |
| `--distance_ip`         | `200`      | `60`             | Internal priming sites closer than this to the isoform 3' end are ignored.                    |
| `--gene_fraction`       | `98%`      | `98%`            | Genes used to learn the read distance distribution (top expressed excluded).                  |
| `--binsize`             | `20`       | `20`             | Bin size (nt) of the read distance distribution.                                              |
| `--em_tol`              | `1e-6`     | `0.01`           | EM convergence tolerance.                                                                     |
| `--em_max_iter`         | `1000`     | `30`             | Maximum number of EM iterations.                                                              |
| `--salmon_decoys`       | `false`    | `false`          | Add the genome as salmon decoys.                                                              |
| `--libtype`             | `A`        | `A`              | Salmon library type.                                                                          |

The defaults that differ from the original SCALPEL are evaluated in the [usage documentation](docs/usage.md) (simulations with ground truth and GSE104556 data). The original settings are reproduced with `--distance_3end 50 --distance_profile 600 --min_mapq 0 --distance_ip 60 --em_tol 0.01 --em_max_iter 30`.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](https://nf-co.re/scalpel/parameters).

## Pipeline output

| Directory          | Description                                                                                          |
| ------------------ | ---------------------------------------------------------------------------------------------------- |
| `dge/`             | **Isoform count matrices** per sample (`<sample>_iDGE/`, sparse 10x format: `matrix.mtx.gz`, `features.tsv.gz` as `gene***isoform`, `barcodes.tsv.gz`). |
| `differential/`    | **Differential isoform usage** (`DIU_table.csv`) and Seurat object of all samples (`iDGE_seurat.RDS`). |
| `filtered_bams/`   | BAM files of the reads used for the quantification (e.g. for coverage tracks).                      |
| `isoform/`         | Selected isoforms per chromosome (`*.annots.gz`, reusable with `--isoform_annotation`) and collapsing summary. |
| `ip_reference/`    | Internal priming reference built from the genome (reusable with `--ip_reference`).                  |
| `probability/`     | Read distance distribution to the 3' end per sample (`.prob`, `.pdf`).                              |
| `salmon/`          | Salmon pseudobulk quantification per sample.                                                        |
| `transcript/`      | Annotation per chromosome with the averaged TPM fractions.                                          |
| `bed/`, `reads/`, `ip/` | Reads per chromosome (filtered, deduplicated), reads mapped to isoforms, after internal priming filtering. |
| `fragment/`, `cell/`, `em/` | Fragment probabilities, cell chunks and EM isoform abundances.                             |
| `pipeline_info/`   | Execution reports, software versions, collapsed isoforms log.                                       |

For more details, please refer to the [output documentation](docs/output.md).

### Downstream analysis

For tutorials on working with SCALPEL output:

- [10X Chromium scRNA-seq example](https://raw.githack.com/plasslab/SCALPEL/refs/heads/dev/docs/chromium.html)
- [Drop-seq scRNA-seq example](https://raw.githack.com/plasslab/SCALPEL/refs/heads/dev/docs/dropseq.html)
- [10X Chromium downstream analysis](https://raw.githack.com/plasslab/SCALPEL/refs/heads/dev/docs/chromium_downstream.html)
- [Drop-seq downstream analysis](https://raw.githack.com/plasslab/SCALPEL/refs/heads/dev/docs/dropseq_downstream.html)

### Testing

`test.qmd` runs the whole test suite (end-to-end runs, regression tests of the plasslab/SCALPEL issues, unit tests of the filtering, collapsing and EM rules) and renders a report:

```bash
quarto render test.qmd
```

## Credits

nf-core/scalpel was originally written by Franz AKE.

We thank the following people for their extensive assistance in the development of this pipeline:

- Marcel Schilling
- Sandra M. Fernández-Moya
- Akshay Jaya Ganesh
- Ana Gutiérrez-Franco
- Lei Li
- Mireya Plass

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#scalpel` channel](https://nfcore.slack.com/channels/scalpel) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/scalpel for your analysis, please cite:

> Franz Ake, Marcel Schilling, Sandra M. Fernández-Moya, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass.
> **Quantification of transcript isoforms at the single-cell level using SCALPEL.**
> _Nat Commun_ 16, 6402 (2025). doi: [10.1038/s41467-025-61118-0](https://doi.org/10.1038/s41467-025-61118-0)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
