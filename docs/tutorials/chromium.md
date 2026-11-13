# SCALPEL application on 10x scRNA‐seq: mouse spermatogenesis



## Introduction

<div style="text-align: center;">

<img src="../vignettes/img/10x_Genomics_logo.svg.png" alt="10x Genomics Logo" width="100px"/> <img src="../vignettes/img/nextflow_logo.png" alt="Nextflow Logo" width="100px"/>

</div>

This vignette demonstrates how to run the nf-core/scalpel pipeline on 10X Genomics single-cell RNA-seq data derived from mouse spermatogenesis, originally published by [Winterpacht A, Lukassen](https://pubmed.ncbi.nlm.nih.gov/30204153/).

The full dataset is publicly available on GEO at accession [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556).

SCALPEL requires as input for each sample the FASTQ files and the Cell Ranger output directory (indexed BAM file and filtered count matrix; Cell Ranger ARC and STARsolo outputs are also supported). For the purpose of this demonstration, we provide a subset of the original data: [reduced_GSE104556](https://zenodo.org/records/17176865/files/reduced_GSE104556.tar.gz?download=1).

``` sh
> wget -O reduced_GSE104556.tar.gz https://zenodo.org/records/17176865/files/reduced_GSE104556.tar.gz?download=1
> tar -xvf reduced_GSE104556.tar.gz
> tree -L 2 reduced_GSE104556
reduced_GSE104556
├── cellranger_samples
│   ├── SRR6129050
│   └── SRR6129051
├── fastq_files
│   ├── SRR6129050_1000000_S1_L001_R1_001.fastq.gz
│   ├── SRR6129050_1000000_S1_L001_R2_001.fastq.gz
│   ├── SRR6129051_1000000_S1_L001_R1_001.fastq.gz
│   └── SRR6129051_1000000_S1_L001_R2_001.fastq.gz
├── SRR6129050_barcodes.txt
└── SRR6129051_barcodes.txt

5 directories, 6 files
```

SCALPEL requires a **GTF genome annotation** and a **FASTA transcript annotation** of your organism, and the **genome FASTA** (used to build the internal priming reference). For the current analysis, we use a reduced annotation reference from GENCODE for Mouse (GENCODE.vM21).

- [reduced_GENCODE.vM21.annotation.gtf](https://zenodo.org/records/17176838/files/reduced_GENCODE.vM21.annotation.gtf.gz?download=1)
- [GENCODE.vM21.transcript.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.transcripts.fa.gz)
- [mm10 genome (UCSC)](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)

``` sh
> wget -O reduced_GENCODE.vM21.annotation.gtf.gz https://zenodo.org/records/17176838/files/reduced_GENCODE.vM21.annotation.gtf.gz?download=1
> wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.transcripts.fa.gz
> wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
```

Complete GTF genome annotation and FASTA transcripts annotation can be downloaded from [GENCODE](https://www.gencodegenes.org) or [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html).

## SCALPEL running

### Samplesheet

To run SCALPEL, you must provide a **samplesheet.csv** file (`--samplesheet`), with a header and one row per sample:

1.  **sample**: sample name (no spaces).
2.  **fastq_1**: R1 FASTQ file.
3.  **fastq_2**: R2 FASTQ file.
4.  **cranger_path**: Cell Ranger output directory of the sample (the BAM file, its index and the filtered barcodes are found automatically).

Paths can be absolute or relative to the samplesheet directory. A sample sequenced in several runs is given as several rows with the same sample name.

``` sh
> cat samplesheet.csv
sample,fastq_1,fastq_2,cranger_path
SRR6129050,reduced_GSE104556/fastq_files/SRR6129050_1000000_S1_L001_R1_001.fastq.gz,reduced_GSE104556/fastq_files/SRR6129050_1000000_S1_L001_R2_001.fastq.gz,reduced_GSE104556/cellranger_samples/SRR6129050
SRR6129051,reduced_GSE104556/fastq_files/SRR6129051_1000000_S1_L001_R1_001.fastq.gz,reduced_GSE104556/fastq_files/SRR6129051_1000000_S1_L001_R2_001.fastq.gz,reduced_GSE104556/cellranger_samples/SRR6129051
```

### Execution

To execute SCALPEL, use the following command:

``` sh
nextflow run plasslab/nf-core-scalpel -r dev \
    -profile conda \
    --samplesheet samplesheet.csv \
    --transcriptome gencode.vM21.transcripts.fa.gz \
    --gtf reduced_GENCODE.vM21.annotation.gtf.gz \
    --genome mm10.fa.gz \
    --outdir results
```

By default, the **internal priming reference is built by the pipeline** from the genome (`--genome mm10.fa.gz`), on the chromosomes of the GTF, in parallel with the salmon quantification. It is saved to `results/ip_reference/ip_reference.track.gz` and can be reused in later runs with `--ip_reference` (see [Pre-built internal priming reference](#pre-built-internal-priming-reference)).

### Options

#### Cell barcodes

You can restrict the analysis to selected cells with a **barcodes_whitelist.csv** file (`--barcodes_whitelist`) giving the barcode file of each sample (one barcode per line, as in the BAM `CB` tag). For example, with the 100 barcodes per sample provided in the `reduced_GSE104556` directory:

``` sh
> cat barcodes_whitelist.csv
sample,barcodes
SRR6129050,reduced_GSE104556/SRR6129050_barcodes.txt
SRR6129051,reduced_GSE104556/SRR6129051_barcodes.txt
```

and add `--barcodes_whitelist barcodes_whitelist.csv` to the command.

#### Cell clusters

To compare cell clusters instead of samples in the differential isoform usage, provide a `--clusters` file (`sample,barcode,cluster`).

#### Pre-built internal priming reference

The internal priming reference built from `--genome` is saved to `<outdir>/ip_reference/`. It can be reused, or replaced by the provided [mm10 reference](https://zenodo.org/records/15664563/files/mm10_polya.track.tar.gz?download=1), with `--ip_reference` (`--genome` is then not needed):

``` sh
    --ip_reference mm10_polya.track.tar.gz
```

The generation from the genome works for any assembly or species (e.g. mm39, T2T).

#### Profiles and resources

Use `-profile docker`, `singularity` or `apptainer` instead of `conda` (add `emulate_amd64` on Apple Silicon with Docker), and `-resume` to restart a run from its completed steps.

> [!TIP]
>
> The salmon steps request 36 GB of memory by default. On a laptop or a small server, cap the resources with a config file given with `-c local.config`:
>
> ``` sh
> > cat local.config
> process {
>     resourceLimits = [cpus: 8, memory: '14.GB', time: '6.h']
> }
> ```

> [!WARNING]
>
> The BAM files, the GTF annotation and the internal priming reference must use the same chromosome names (e.g. ‘chr1’, …, ‘chrM’ vs ‘1’, …, ‘MT’). The pipeline checks this before processing the reads and stops with an explicit message in case of mismatch. Use references from the same source as the BAM files.

### Main parameters

| Parameter | Default | Description |
|----|----|----|
| `--distance_3end` | `100` | Maximum distance (nt) between the 3’ ends of collapsed isoforms. |
| `--distance_profile` | `500` | Window (nt from the 3’ end) where collapsed isoforms have identical exons. |
| `--distance_reads` | `600` | Maximum distance (nt) of a read from the isoform 3’ end. |
| `--min_mapq` | `255` | Minimum mapping quality (255: uniquely mapped reads). |
| `--distance_ip` | `200` | Internal priming sites closer to the 3’ end are ignored. |
| `--gene_fraction` | `98%` | Genes used to learn the read distance distribution. |
| `--binsize` | `20` | Bin size (nt) of the read distance distribution. |

See the [usage documentation](https://github.com/plasslab/nf-core-scalpel/blob/dev/docs/usage.md) for all options.

### Results

``` sh
> nextflow run plasslab/nf-core-scalpel -r dev -profile conda -c local.config --samplesheet samplesheet.csv --transcriptome gencode.vM21.transcripts.fa.gz --gtf reduced_GENCODE.vM21.annotation.gtf.gz --genome mm10.fa.gz --outdir results

------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/scalpel 1.0.0dev
------------------------------------------------------
Input/output options
  samplesheet        : samplesheet.csv
  gtf                : reduced_GENCODE.vM21.annotation.gtf.gz
  genome             : mm10.fa.gz
  transcriptome      : gencode.vM21.transcripts.fa.gz
  outdir             : results
!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x
* Software dependencies
    https://github.com/nf-core/scalpel/blob/master/CITATIONS.md
Pipeline options
  sequencing_platform : chromium
  libtype             : A
  distance_profile    : 500
  distance_reads      : 600
  min_mapq            : 255
  distance_3end       : 100
  distance_ip         : 200
  ip_reference        : generated from --genome
  salmon_decoys       : false
  gene_fraction       : 98%
  binsize             : 20
  em_tol / max_iter   : 0.000001 / 1000
  barcodes_whitelist  : not provided
------------------------------------------------------
[PIPELINE] nf-core/scalpel 1.0.0dev | profile=conda
[PROCESS a2/ee760b] READS_PROCESSING:IP_REFERENCE_GENERATION (mm10.fa.gz)
[PROCESS b0/d506a1] ANNOTATION_PROCESSING:SALMON_INDEX (<(gunzip -c gencode.vM21.transcripts.fa.gz))
[PROCESS 4c/317ccd] ANNOTATION_PROCESSING:SALMON_QUANT (SRR6129050)
[PROCESS a1/59ca25] ANNOTATION_PROCESSING:SALMON_QUANT (SRR6129051)
[PROCESS 29/8271ed] ANNOTATION_PROCESSING:TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING (reduced_GENCODE.vM21.annotation.gtf.gz)
[PROCESS 17/9de559] ANNOTATION_PROCESSING:ISOFORM_SELECTION (chr18)
...
[PROCESS 95/4d4884] ISOFORM_QUANTIFICATION:DGE_GENERATION (SRR6129050)
[PROCESS ac/423c15] ISOFORM_QUANTIFICATION:DGE_GENERATION (SRR6129051)
[PROCESS de/7e4ebb] APA_CHARACTERIZATION:DIFFERENTIAL_ISOFORM_USAGE (2 samples)
-[nf-core/scalpel] Pipeline completed successfully-
[SUCCESS] completed=46 failed=0 cached=0
```

Following the successful execution of SCALPEL, the output directory (`results`) contains:

- **differential/iDGE_seurat.RDS**:  
  The Seurat object containing the isoform expression data (iDGE) of all samples (one counts layer per sample, cells named `<sample>_<barcode>`). This object can be used for downstream analysis and visualization in R.

- **differential/DIU_table.csv**:  
  A table summarizing the differential isoform usage (DIU) events between the samples (or the clusters given with `--clusters`).

- **dge/sampleXXX_iDGE/**:  
  The isoform count matrix (iDGE) of each sample, in sparse 10x format (`matrix.mtx.gz`, `features.tsv.gz` with `gene***isoform` features, `barcodes.tsv.gz`), readable with `Seurat::Read10X()`.

- **filtered_bams/sampleXXX_filtered.bam(.bai)**:  
  The reads used for the quantification (duplicated reads and internal priming associated reads discarded, only reads passing the SCALPEL filters).

- **isoform/**:  
  The selected isoforms, reusable with `--isoform_annotation` to quantify new samples on the same isoforms.

- **ip_reference/**:  
  The internal priming reference built from the genome (reusable with `--ip_reference`).

- **probability/**:  
  The read distance distribution to the 3’ end of each sample.

- **pipeline_info/**:  
  Execution reports and software versions.

The other directories contain intermediate files of each step (see the [output documentation](https://github.com/plasslab/nf-core-scalpel/blob/dev/docs/output.md)).

``` sh
> tree -L 2 results
results
├── bed
│   ├── SRR6129050
│   └── SRR6129051
├── cell
│   ├── SRR6129050_cell_chunk_1.fst
│   ├── SRR6129050_cell_chunk_2.fst
│   ├── SRR6129050_cell_chunk_3.fst
│   ├── SRR6129051_cell_chunk_1.fst
│   ├── SRR6129051_cell_chunk_2.fst
│   └── SRR6129051_cell_chunk_3.fst
├── dge
│   ├── SRR6129050_iDGE
│   └── SRR6129051_iDGE
├── differential
│   ├── DIU_table.csv
│   └── iDGE_seurat.RDS
├── em
│   ├── SRR6129050_cell_chunk_1_em_results.fst
│   ├── SRR6129050_cell_chunk_2_em_results.fst
│   ├── SRR6129050_cell_chunk_3_em_results.fst
│   ├── SRR6129051_cell_chunk_1_em_results.fst
│   ├── SRR6129051_cell_chunk_2_em_results.fst
│   └── SRR6129051_cell_chunk_3_em_results.fst
├── filtered_bams
│   ├── SRR6129050_filtered.bam
│   ├── SRR6129050_filtered.bam.bai
│   ├── SRR6129051_filtered.bam
│   └── SRR6129051_filtered.bam.bai
├── fragment
│   ├── SRR6129050_SRR6129050_chr10.ip_filtered.frag.fst
│   ├── SRR6129050_SRR6129050_chr17.ip_filtered.frag.fst
│   ├── SRR6129051_SRR6129051_chr10.ip_filtered.frag.fst
│   └── SRR6129051_SRR6129051_chr17.ip_filtered.frag.fst
├── ip
│   ├── SRR6129050_chr10.ip_filtered.fst
│   ├── SRR6129050_chr17.ip_filtered.fst
│   ├── SRR6129051_chr10.ip_filtered.fst
│   └── SRR6129051_chr17.ip_filtered.fst
├── ip_reference
│   └── ip_reference.track.gz
├── isoform
│   ├── chr10.annots.gz
│   ├── chr12.annots.gz
│   ├── chr17.annots.gz
│   ├── chr18.annots.gz
│   └── stats.txt
├── pipeline_info
│   ├── collapsed_isoforms.log
│   ├── execution_report_2026-11-12_01-16-10.html
│   ├── execution_timeline_2026-11-12_01-16-10.html
│   ├── execution_trace_2026-11-12_01-16-10.txt
│   ├── nf_core_scalpel_software_versions.yml
│   ├── params_2026-11-12_01-16-12.json
│   └── pipeline_dag_2026-11-12_01-16-10.html
├── probability
│   ├── SRR6129050.pdf
│   ├── SRR6129050.prob
│   ├── SRR6129051.pdf
│   └── SRR6129051.prob
├── reads
│   ├── SRR6129050_chr10.fst
│   ├── SRR6129050_chr17.fst
│   ├── SRR6129051_chr10.fst
│   └── SRR6129051_chr17.fst
├── salmon
│   ├── salmon
│   ├── SRR6129050
│   ├── SRR6129050_lib_format_counts.json
│   ├── SRR6129050_meta_info.json
│   ├── SRR6129051
│   ├── SRR6129051_lib_format_counts.json
│   └── SRR6129051_meta_info.json
└── transcript
    ├── chr10.tsv
    ├── chr12.tsv
    ├── chr17.tsv
    └── chr18.tsv
```

## Downstream analysis

For more information, refer to [SCALPEL_downstream_analysis](https://raw.githack.com/plasslab/nf-core-scalpel/refs/heads/dev/docs/vignettes/chromium_downstream.html) for downstream analysis of the SCALPEL output files.

email: <a href="mailto:aerod7710@gmail.com"> fake@idibell.cat </a>
