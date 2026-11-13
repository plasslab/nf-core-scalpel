# Installation

## Pipeline

SCALPEL runs with [Nextflow](https://www.nextflow.io/) (>= 26.04) and a software profile providing the tools of each step:

| Profile                   | Use                                                                          |
| ------------------------- | ---------------------------------------------------------------------------- |
| `docker`                  | Local computer (add `emulate_amd64` on Apple Silicon).                       |
| `singularity`, `apptainer`| HPC clusters and servers.                                                    |
| `conda`                   | Without containers (slower first run: environments are created).             |

If you are new to Nextflow and nf-core, see the [nf-core installation page](https://nf-co.re/docs/usage/installation). The pipeline is
downloaded by Nextflow at the first run:

```bash
nextflow run plasslab/nf-core-scalpel -r dev -profile test,docker --outdir results_test
```

The `test` (10x Chromium) and `test_dropseq` profiles run the pipeline on small public datasets in a few minutes.

### Resources

The default resources of each step follow the nf-core labels (e.g. salmon asks 36 GB of memory). On a computer with less memory or CPUs,
cap the resources with a configuration file passed with `-c`:

```groovy title="local.config"
process {
    resourceLimits = [cpus: 8, memory: 16.GB, time: 24.h]
}
```

```bash
nextflow run plasslab/nf-core-scalpel -r dev -profile docker -c local.config ...
```

## scalpelR

The downstream analysis functions (differential isoform usage between groups of cells, isoform usage plots, read coverage, 3' UTR length
per cell, 3' UTR shift heatmaps) are provided by the R package [scalpelR](https://github.com/Franzx7/scalpelR):

```r
install.packages("remotes")
remotes::install_github("Franzx7/scalpelR@dev")
```

scalpelR uses Seurat (>= 5), data.table, Matrix, ggplot2 and patchwork; `CoveragePlot` also needs
[samtools](http://www.htslib.org/) in the `PATH`. See the [scalpelR reference](scalpelR/index.md).
