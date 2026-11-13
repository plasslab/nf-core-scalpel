# FAQ

## Inputs

**Which data can SCALPEL analyse?**
3'-tagged scRNA-seq: 10x Genomics Chromium 3' (Cell Ranger, Cell Ranger ARC or STARsolo outputs) and Drop-seq (Drop-seq tools BAM and
DGE count matrix, `--sequencing_platform dropseq`). The FASTQ files are also needed for the salmon pseudobulk quantification
([samplesheet](usage.md#samplesheet-input)).

**The pipeline stops because the BAM and the GTF share no chromosome name.**
The BAM files, the GTF and the genome must use the same chromosome names (e.g. GENCODE/UCSC `chr1` vs Ensembl `1`). Use references from
the same source as the reference of your BAM files ([details](usage.md#chromosome-names-must-match)).

**Can I add samples to a previous run?**
Yes: `--isoform_annotation <outdir>/isoform` reuses the isoforms selected in the previous run, so that the new samples are quantified on
the same isoforms ([details](usage.md#adding-samples-later-isoform-annotation)).

**Should I restrict the analysis to selected cells?**
Recent Cell Ranger versions call many low-RNA barcodes. Filter the cells at the gene level (number of genes per cell, mitochondrial
reads) and either give the kept barcodes with `--barcodes_whitelist`, or keep them in the downstream analysis (see the
[10x downstream tutorial](tutorials/chromium_downstream.md)).

## Differential isoform usage

**Samples or cell clusters?**
Without `--clusters`, the samples are compared. With one library per condition, biological and technical (batch) differences cannot be
separated. Compare cell clusters or cell types with `--clusters` (a `sample`, `barcode`, `cluster` file): with clusters shared by several
samples, the permutation test is stratified by sample and the replication of the changes is assessed
([details](usage.md#cell-clusters-clusters-for-differential-isoform-usage)).

**`chisq` or `permutation` test?**
The chi-squared test on pooled counts is fast, but ignores the variability between cells and calls many changes on large datasets. The
permutation test (`--diu_test permutation`) uses the cells as units ([details](usage.md#differential-isoform-usage-test-diu-test)).

**What does `replicated` mean in `DIU_table.csv`?**
When the compared clusters are shared by several samples, the change of the most changing isoform of each gene is measured in each sample
with enough UMIs. A gene is significant only if its change is `replicated` (same direction in each evaluable sample, at least 2 such
samples). `NA`: replication not applicable (samples compared, or each cluster in its own sample).

**What are `tandem_3UTR` and `alternative_terminal_exon` events?**
`tandem_3UTR`: the tested isoforms share their terminal exon (alternative polyadenylation in the last exon), with a `direction`
(`shortening` / `lengthening` of the 3' UTRs of group1 vs group2). `alternative_terminal_exon`: the isoforms end in different exons; no
3' UTR direction is given, `switch` tells whether the dominant isoform differs between the groups.

## Resources

**A step fails with an out-of-memory error.**
Cap or raise the resources with a configuration file passed with `-c` (see [Installation](installation.md#resources)) and restart the
run with `-resume`: the completed steps are not run again.

**How long does a run take?**
On the full datasets of the tutorials (2 samples, 64 CPUs): about 1 h 10 min for 10x Chromium (GSE104556) and 40 min for Drop-seq
(GSE268222). Re-running only the differential isoform usage with `-resume` (e.g. new `--clusters`) takes a few minutes.

## Downstream analysis

**How do I read the isoform-based Seurat object?**
`differential/iDGE_seurat.RDS` has one counts layer per sample (Seurat v5) and cells named `<sample>_<barcode>`: join the layers with
`SeuratObject::JoinLayers()` before the analysis. Features are named `gene***isoform`.

**Is the 3' UTR length per cell affected by the sequencing depth?**
`scalpelR::weightedUTR_cell` compares each gene with its own weighted 3' UTR length in all the cells (`method = "centered"`, default),
so the value does not depend on which genes are detected in a cell. The former metric (`method = "length"`, mean 3' UTR length of the
genes detected in the cell) increases with the sequencing depth.
