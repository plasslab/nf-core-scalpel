# SCALPEL: isoforms in single cells

**SCALPEL** quantifies transcript isoforms and alternative polyadenylation (APA) at the single-cell level from 3'-tagged scRNA-seq
data (10x Genomics Chromium 3', Drop-seq). It is distributed as a Nextflow pipeline (**nf-core/scalpel**) and an R package for the
downstream analysis (**scalpelR**).

```{image} images/SCALPEL_diagram.png
:alt: SCALPEL pipeline overview
:width: 800px
:align: center
```

1. **Annotation processing**: salmon pseudobulk quantification and selection of the isoforms distinguishable by 3'-tagged reads
   (isoforms with close 3' ends and identical exons near the 3' end are collapsed).
2. **Reads processing**: reads of the filtered cells mapped to the selected isoforms, internal priming artifacts removed.
3. **Isoform quantification**: EM estimation of the isoform abundances per cell and gene, isoform-by-cell count matrices (iDGE).
4. **APA characterization**: differential isoform usage between samples or cell clusters.

## Getting started

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`download;1.5em` Installation
:link: installation
:link-type: doc

Nextflow, container profiles and the scalpelR package.
:::

:::{grid-item-card} {octicon}`play;1.5em` Usage
:link: usage
:link-type: doc

Samplesheet, references, cell clusters, differential isoform usage and parameters.
:::

:::{grid-item-card} {octicon}`file-directory;1.5em` Output
:link: output
:link-type: doc

Isoform count matrices, Seurat object, differential isoform usage table.
:::

:::{grid-item-card} {octicon}`mortar-board;1.5em` Tutorials
:link: tutorials/index
:link-type: doc

Full runs and downstream analyses on 10x Chromium and Drop-seq data.
:::
::::

## Quick test

```bash
nextflow run plasslab/nf-core-scalpel -r dev -profile test,docker --outdir results_test
```

Use `conda`, `singularity` or `apptainer` instead of `docker`, `test_dropseq` for the Drop-seq test data, and add `emulate_amd64` on
Apple Silicon with Docker.

## Citation

If you use SCALPEL for your analysis, please cite:

> Franz Ake, Marcel Schilling, Sandra M. Fernández-Moya, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass.
> **Quantification of transcript isoforms at the single-cell level using SCALPEL.**
> _Nat Commun_ 16, 6402 (2025). doi: [10.1038/s41467-025-61118-0](https://doi.org/10.1038/s41467-025-61118-0)

```{toctree}
:hidden:
:maxdepth: 2

installation
usage
output
tutorials/index
scalpelR/index
faq
changelog
```
