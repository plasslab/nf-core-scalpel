#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(Matrix, quietly = TRUE))


parser = ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('--dge', type="character", help='path of DGE directory (10x filtered_feature_bc_matrix)')
parser$add_argument('--pred', type="character", help='path of predicted isoforms')
parser$add_argument('--output', type="character", help="path of output directory (10x format iDGE)")
args = parser$parse_args()


# File opening
# ++++++++++++

#merge all em files
list_em_files = list.files(args$pred, full.names = TRUE, pattern = "*em_results.fst")

#load cell isoform probabilities
message("Loading isoform predictions...")
preds = lapply(list_em_files, function(x) {
    tmp = read_fst(x, as.data.table = TRUE)
    return(tmp)
}) %>% rbindlist()

#load DGE count matrix (sparse, no dense conversion #38)
message("Loading DGE matrix...")
dge_counts = Matrix::readMM(file.path(args$dge, "matrix.mtx.gz")) %>% as("TsparseMatrix")
dge_genes = fread(file.path(args$dge, "features.tsv.gz"), header = FALSE)$V2
dge_bcs = fread(file.path(args$dge, "barcodes.tsv.gz"), header = FALSE)$V1
dge = data.table(
    gene_idx = dge_counts@i + 1L,
    bc_idx = dge_counts@j + 1L,
    counts = dge_counts@x
)
rm(dge_counts)


# Merging DGE & isoform predictions
# +++++++++++++++++++++++++++++++++
message("Merging DGE & isoform predictions...")

#- first gene symbol occurrence is kept (as with Read10X rownames)...
preds[, gene_idx := match(gene_name, dge_genes)]
preds[, bc_idx := match(bc, dge_bcs)]
preds = dge[preds, on = .(gene_idx, bc_idx)]
preds[is.na(counts), counts := 0]

#- isoform counts...
preds[, tr_counts := as.numeric(rel_abund * counts)]
preds[, gene_transcript := paste0(gene_name, "***", transcript_name)]
preds = unique(preds, by = c("gene_transcript", "bc"))


# Sparse iDGE matrix (isoforms x cells)
# +++++++++++++++++++++++++++++++++++++
message("Building sparse iDGE matrix...")
features = sort(unique(preds$gene_transcript))
barcodes = sort(unique(preds$bc))
preds = preds[tr_counts != 0]

idge = Matrix::sparseMatrix(
    i = match(preds$gene_transcript, features),
    j = match(preds$bc, barcodes),
    x = preds$tr_counts,
    dims = c(length(features), length(barcodes))
)
message("iDGE: ", length(features), " isoforms x ", length(barcodes), " cells, ", length(idge@x), " non-zero entries")


#writing (10x format, readable with Seurat::Read10X(gene.column = 1))
dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
Matrix::writeMM(idge, file.path(args$output, "matrix.mtx"))
R.utils::gzip(file.path(args$output, "matrix.mtx"), overwrite = TRUE)
fwrite(data.table(features), file.path(args$output, "features.tsv.gz"), col.names = FALSE)
fwrite(data.table(barcodes), file.path(args$output, "barcodes.tsv.gz"), col.names = FALSE)
