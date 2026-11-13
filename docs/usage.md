# nf-core/scalpel: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/scalpel/usage](https://nf-co.re/scalpel/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

SCALPEL quantifies isoforms (alternative 3' ends) from 3'-tagged scRNA-seq data
(10x Chromium 3', Drop-seq). It needs, per sample, the raw FASTQ files (pseudobulk
isoform quantification with salmon) and the aligned, cell-barcoded BAM file with its
filtered count matrix (Cell Ranger, Cell Ranger ARC or STARsolo output).

## Samplesheet input

```bash
--samplesheet '[path to samplesheet.csv]'
```

A comma-separated file with a header and 4 columns:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cranger_path
SRR6129050,reads/SRR6129050_R1.fastq.gz,reads/SRR6129050_R2.fastq.gz,cellranger/SRR6129050
SRR6129051,reads/SRR6129051_R1.fastq.gz,reads/SRR6129051_R2.fastq.gz,cellranger/SRR6129051
```

| Column         | Description                                                                                                                  |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Sample name (no spaces). Used in all outputs and in the `--clusters` file.                                                   |
| `fastq_1`      | R1 FASTQ (`.fastq.gz`/`.fq.gz`). For 10x, R1 holds barcode + UMI.                                                             |
| `fastq_2`      | R2 FASTQ (cDNA read).                                                                                                        |
| `cranger_path` | 10x: the sample output **directory** (see layouts below). Drop-seq (`--sequencing_platform dropseq`): the tagged BAM **file**. |
| `dge`          | Drop-seq only (required): the Drop-seq tools DGE count matrix (`.txt`/`.txt.gz`). Leave out for 10x.                         |

- Paths can be absolute or **relative to the samplesheet directory** (not the launch directory).
- Every path is checked when the pipeline starts; a missing file stops the run with a
  message naming the sample and the column.
- Do not use `NA` columns: the FASTQ files are required (pseudobulk quantification).

### Several runs (FASTQ pairs) of one sample

A sample sequenced in several runs or lanes is given as several rows with the same `sample` name,
one FASTQ pair per row, all pointing to the same `cranger_path` (and `dge` for Drop-seq): all its
FASTQ pairs are quantified together by salmon.

```csv
sample,fastq_1,fastq_2,cranger_path
S1,S1_L001_R1.fastq.gz,S1_L001_R2.fastq.gz,cellranger/S1
S1,S1_L002_R1.fastq.gz,S1_L002_R2.fastq.gz,cellranger/S1
```

### Supported sample layouts (`cranger_path`)

The pipeline detects the layout of each sample directory:

| Layout                      | BAM                                  | Count matrix                         |
| --------------------------- | ------------------------------------ | ------------------------------------ |
| Cell Ranger `count`         | `outs/possorted_genome_bam.bam`      | `outs/filtered_feature_bc_matrix/`   |
| Cell Ranger ARC (multiome)  | `outs/gex_possorted_bam.bam`         | `outs/filtered_feature_bc_matrix/`   |
| STARsolo                    | `Aligned.sortedByCoord.out.bam`      | `Solo.out/Gene/filtered/`            |

- **Multiome (GEX + ATAC)**: only the `Gene Expression` features of the count matrix are
  used, the ATAC `Peaks` features are discarded automatically. No renaming of files is needed.
- **STARsolo**: run STAR with `--outSAMtype BAM SortedByCoordinate --outSAMattributes CB UB`
  (cell barcode / UMI tags are required). A missing BAM index is created automatically.
- **Drop-seq** (`--sequencing_platform dropseq`): see below.

### Drop-seq samplesheet

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,cranger_path,dge
S1,reads/S1_R1.fastq.gz,reads/S1_R2.fastq.gz,bams/S1.bam,dges/S1.dge.txt.gz
```

- `cranger_path`: coordinate-sorted BAM from Drop-seq tools with the `XC` (cell barcode) and
  `XM` (UMI) tags; an index is created if missing.
- `dge`: the Drop-seq tools digital expression matrix (`DigitalExpression`): tab-separated,
  first column `GENE` (gene symbols), then one column per cell barcode.
- The cells analysed are the barcodes of the DGE header (as the filtered barcodes for 10x),
  unless a `--barcodes_whitelist` is given.

## Reference files

| Parameter         | Description                                                                                                  |
| ----------------- | ------------------------------------------------------------------------------------------------------------ |
| `--transcriptome` | Transcript FASTA (e.g. `gencode.vM21.transcripts.fa.gz`). **Not** the genome FASTA.                           |
| `--gtf`           | Matching annotation GTF (same release as the transcriptome).                                                  |
| `--genome`        | Genome FASTA, used to build the internal priming reference (default, see below).                             |
| `--ip_reference`  | Optional. Pre-built internal priming sites, used instead of building them from `--genome` (see below).        |
| `--salmon_decoys` | Optional. Adds the `--genome` sequences as salmon decoys. Not used by the original SCALPEL (default: false).  |

### Chromosome names must match

The BAM files (i.e. the reference used by Cell Ranger/STARsolo), the GTF and the internal
priming reference must use the **same chromosome names** (e.g. GENCODE/UCSC `chr1` vs
Ensembl `1`). The pipeline checks this before processing the reads and stops with an
explicit message if the BAM and the GTF share no chromosome name (instead of silently
skipping the downstream steps). Use references from the same source as your BAM files;
renaming by hand (adding/removing `chr`, `M`/`MT`) is error-prone.

### Internal priming reference

The internal priming sites are provided in two equivalent formats depending on the
download: `<genome>.polyA.track` (e.g. `mm10_polya.track.gz`) or
`<genome>_ipriming_sites.bed(.tar.gz)` (e.g. `hg38_ipriming_sites.bed.tar.gz`). Both contain
the same information (tab-separated: chromosome, start, end, ., ., strand) and can be given
to `--ip_reference` as plain, `.gz` or `.tar.gz` files. The reference is split by chromosome
once at the start of the run.

**Generation from the genome (default).** Without `--ip_reference`, the internal priming sites
are detected on the annotated chromosomes of `--genome` (plain or `.gz` FASTA, same chromosome
names as the GTF), with the rule of the SCALPEL mm10/hg38 references: genomic windows of 10 nt
with at least 7 A, or starting with a stretch of 6 A, overlapping windows merged and trimmed to
their first A (`+` strand; `-` strand: same rule on T, reverse complement). On mouse chr17, the
generated sites match the SCALPEL `mm10_polya.track` reference at 99.8% (identical sites) and
99.99% (base-level Jaccard index). The reference is saved to
`<outdir>/ip_reference/ip_reference.track.gz` and can be reused with `--ip_reference`.
It starts with the run (chromosomes read from the GTF) and runs in parallel with salmon: about
4 min for a mouse genome (2.7 Gb, 4-6 CPUs, ~1.3 GB memory per 100 Mb chromosome).

```bash
nextflow run nf-core/scalpel \
   -profile <docker/singularity/.../institute> \
   --samplesheet samplesheet.csv \
   --transcriptome /path/to/transcriptome.fa \
   --gtf /path/to/annotation.gtf \
   --genome /path/to/genome.fa \
   --outdir <OUTDIR>
```

**Filtering rule (original SCALPEL).** For each isoform, the internal priming site closest to
its 3' end, at least `--distance_ip` nt away (default 200), is located on its exons. Every
fragment (cell + UMI) lying entirely upstream of that site is discarded for this isoform, as it
more likely originates from oligo-dT priming on that A-rich stretch than from the genuine 3'
end. Sites closer than `--distance_ip` to the 3' end are ignored. The rule applies to both
strands.

The default of 200 nt (original SCALPEL: 60) keeps sites inside the bulk of an isoform's own
read distribution (median 110-130 nt from the 3' end) from discarding its genuine reads. On the
GSE104556 test data, 200 nt cuts fragments removed despite a 3' end peak on their isoform by
~78% compared to 60 nt, at the cost of skipping sites closer than 200 nt. Use
`--distance_ip 60` to reproduce the original behaviour.

## Reads used (`--min_mapq`)

By default only uniquely mapped reads are used (`--min_mapq 255`, the MAPQ of unique alignments
in Cell Ranger, STARsolo and STAR / Drop-seq BAMs), i.e. the reads counted in the gene count
matrix. The isoform counts are the EM isoform proportions x the gene counts: multi-mapped reads,
not counted in the gene matrix, can come from another copy of a paralog family and bias the
proportions. On GSE104556, 25% of the reads used on chrX/chrY (15% of the chr17 reads) are
multi-mapped. Use `--min_mapq 0` to keep all reads, as the original SCALPEL.

## Isoform selection (collapsing)

Within each gene, isoforms whose 3' ends are within `--distance_3end` (default 100 nt) of each other
are compared over the last `--distance_profile` nt (default 500) from the furthest 3' end. Isoforms
with identical exons in this window (same genomic splice sites) cannot be told apart by 3'-tagged
reads and are collapsed into the most expressed one (pseudobulk TPM), whose terminal exon is
extended to the furthest 3' end of the group (`isoform/*.annots.gz`, column `transcript_collapseds`).
Exons of the same length at different genomic positions are different isoforms (not collapsed).

The collapsing window (`--distance_profile`, 500 nt) and the maximum read distance to the isoform
3' end (`--distance_reads`, 600 nt) are separate: the original SCALPEL used 600 nt for both. On
GSE104556, a structural difference 500-600 nt upstream of the 3' end is covered by 3-7% of the
reads, so such isoforms were mostly split by their TPM prior (250 isoforms, 1.3%, now collapsed),
whereas reads keep the 3' signal up to ~600 nt (further upstream, the read density is gene body
background). Use `--distance_profile 600` to reproduce the original isoform selection.

The default `--distance_3end` is 100 nt (original SCALPEL: 50 nt). Simulating two isoforms with
identical exons and 3' ends 50-100 nt apart with the GSE104556 read distribution, the per-cell
fraction estimates are no better than an uninformative guess and pooled estimates are off by
6-15 points, so such isoforms are collapsed (276 more isoforms, 1.4%, on GSE104556). Use
`--distance_3end 50 --distance_profile 600` to reproduce the original isoform selection.

## Read distance distribution

The probability of a read given its distance to the isoform 3' end is learned per sample from
the genes with a single selected isoform (excluding the top `--gene_fraction` expressed genes),
in bins of `--binsize` nt. When it is learned from fewer than 100 reads (very small / test
datasets), the pipeline warns (log and `probability/<sample>.warning.txt`): many distances are
then never observed and get a null probability, and the EM relies on the isoform TPM fractions.

## Read distance distribution

The probability of a read given its distance to the isoform 3' end is learned per sample from
the genes with a single selected isoform (excluding the top `--gene_fraction` expressed genes),
in bins of `--binsize` nt. When it is learned from fewer than 100 reads (very small / test
datasets), the pipeline warns (log and `probability/<sample>.warning.txt`): many distances are
then never observed and get a null probability, and the EM relies on the isoform TPM fractions.

## Isoform quantification (EM)

For each cell and gene, the isoform abundances are estimated by an EM algorithm from the
fragment probabilities (read distances to the isoform 3' ends x isoform TPM fraction). The EM runs
until no abundance changes by more than `--em_tol` (default `1e-6`, at most `--em_max_iter`
1000 iterations). The original SCALPEL stopped at a change < 0.01 (30 iterations), before
isoforms that are never the best explanation of any fragment (e.g. with a 3' end further from all
the reads than another isoform) reach 0: on GSE104556, 3,271 such estimates kept > 5% of the
abundance (15 once converged). Use `--em_tol 0.01 --em_max_iter 30` to reproduce it.

## Optional inputs

### Cell barcodes whitelist (`--barcodes_whitelist`)

By default all cells of the filtered count matrix are used. To restrict the analysis to a
subset of cells (e.g. after your own QC), provide a CSV mapping each sample to a barcode file:

```csv title="barcodes_whitelist.csv"
sample,barcodes
SRR6129050,barcodes/SRR6129050_barcodes.txt
SRR6129051,barcodes/SRR6129051_barcodes.txt
```

Each barcode file has **one barcode per line, no header, no other column**, written as in
the BAM `CB` tag (e.g. `AAACCTGAGCTTATCG-1`). Paths are relative to the CSV. Subsetting the
final `iDGE_seurat.RDS` (e.g. to a cluster) is also valid and does not require re-running.

### Cell clusters (`--clusters`) for differential isoform usage

Without `--clusters`, samples are compared. To compare clusters / cell types / conditions,
provide a tab- or comma-separated file with 3 columns (header optional):

```tsv title="clusters.tsv"
sample	barcode	cluster
SRR6129050	AAACCTGAGCTTATCG-1	Tcells
SRR6129050	AAACCTGGTTGAGTTC-1	Bcells
SRR6129051	AAACCTGTCAACGAAA-1	Tcells
```

- column 1: the sample name **as in the samplesheet**;
- column 2: the cell barcode **as in the BAM / count matrix** (keep the `-1` suffix);
- column 3: the group to compare (at least 2 distinct values).

The pipeline reports how many cells were matched, ignores cells without a cluster and stops
with an explicit message if no cell matches (usually a sample name or barcode suffix mismatch).

### Differential isoform usage test (`--diu_test`)

For each gene, the isoform counts are aggregated per group and the isoform usage (fraction of the
gene counts) is compared between groups: 2 groups are compared with each other, more than 2 groups
each against the rest of the cells. Isoforms below 5% of the gene counts in both groups are not tested
(`thr.abund`); genes are tested when at least 2 isoforms remain and all expected counts are >= 5.

| `--diu_test`            | p-value                                                                                                                              |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `chisq` (default)       | X2 test on the aggregated counts (original SCALPEL). UMIs are treated as independent observations.                                     |
| `permutation`           | Same X2 statistic, null distribution from permutations of the **cells** between the compared groups (`--diu_permutations`, default 10000). |

The X2 test ignores that isoform usage varies from cell to cell: with thousands of cells, small
changes become highly significant (on simulated cells with identical usage in both groups but
cell-to-cell variability, 60% of the genes have p < 0.05 with `chisq`, 4% with `permutation`).
The permutation test uses the cells as units, needs no replicates, and with `--clusters` permutes
the cells **within each sample**, so that a difference between samples (batch) is not reported as
a difference between clusters that are unevenly distributed across samples. Samples compared to each
other (no `--clusters`) are single observations: differences between them can be technical or
biological and cannot be separated without replicates.

The same cell permutations are used for all genes (one sparse matrix product per batch of
permutations), and they stop early for a gene once 20 permuted statistics exceed the observed one, so
only genes with small p-values run the maximum number of permutations (smallest p-value 1/(n+1)).
On GSE104556 (2 samples, ~11,700 cells, 7,000 genes tested), the differential isoform usage step
takes ~1.5 min with `--diu_test permutation` (6 CPUs).

For both tests, the p-values are adjusted (BH) once per gene and comparison, and a gene is
`significant` when its adjusted p-value is < 0.05 **and** the usage of at least one tested isoform
changes by at least `--diu_min_delta` (default 0.1; isoforms with a smaller change: `thr.var`).

```bash
nextflow run nf-core/scalpel ... --clusters clusters.tsv --diu_test permutation
```

### Adding samples later (`--isoform_annotation`)

SCALPEL selects the isoforms once from **all** samples (pseudobulk). To keep the isoform set
fixed when adding samples to an existing analysis, reuse the selected isoforms of the first
run (`<outdir>/isoform/*.annots.gz`):

```bash
nextflow run nf-core/scalpel --samplesheet new_samples.csv --isoform_annotation results_run1/isoform \
    --ip_reference results_run1/ip_reference/ip_reference.track.gz --outdir results_run2 -profile conda
```

Annotation processing (salmon pseudobulk, isoform selection) is skipped (`--transcriptome`,
`--gtf` not needed), and the new samples are quantified on exactly the same isoforms. The FASTQ
files are then not used: the samplesheet columns `fastq_1`/`fastq_2` are optional, so samples can
be given with their BAM / count matrix only (`sample,cranger_path`). Note that
the isoform set then reflects only the samples of the first run; re-running all samples
together remains the only way for the annotation to use the full dataset.

## Resources & large datasets

- The salmon pseudobulk step is usually the longest. Give it enough CPUs/memory
  (`process_high` label) and, on network file systems, prefer a local work directory
  (`-work-dir`). Runs interrupted by out-of-memory/disk errors resume from the last
  completed tasks with `-resume`.
- For 10x data, R1 (barcode + UMI, shorter than the salmon k-mer size) does not map; salmon
  quantifies the R2 reads (reported as orphans). This is expected.
- The isoform count matrices (iDGE) are built and written as sparse matrices, so large
  datasets (> 2^31 isoform x cell entries) are supported without splitting the cells.

## Test profile

```bash
nextflow run . -profile test,conda --outdir results_test
```

Runs on a light dataset (GSE104556, chr17:14-16Mb, 2 samples; `tests/data/light`) in a few
minutes. `-profile test_dropseq,conda` runs the same data in Drop-seq format
(`tests/data/light_dropseq`). `test.qmd` (Quarto) runs it and the regression tests of the reported issues:

```bash
quarto render test.qmd
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scalpel --samplesheet ./samplesheet.csv --transcriptome transcripts.fa.gz --gtf annotation.gtf.gz --genome genome.fa.gz --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scalpel -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scalpel
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scalpel releases page](https://github.com/nf-core/scalpel/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
