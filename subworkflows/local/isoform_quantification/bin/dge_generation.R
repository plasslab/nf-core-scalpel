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
if(dir.exists(args$dge)){
    #- 10x format directory (Cell Ranger: *.gz files, STARsolo: plain files, Cell Ranger v2: genes.tsv)...
    dge_file = function(names){
        candidates = file.path(args$dge, c(paste0(names, ".gz"), names))
        if(!any(file.exists(candidates))) stop("Missing ", paste(names, collapse = "/"), " in ", args$dge)
        return(candidates[file.exists(candidates)][1])
    }
    dge_counts = Matrix::readMM(dge_file("matrix.mtx")) %>% as("TsparseMatrix")
    dge_features = fread(dge_file(c("features.tsv", "genes.tsv")), header = FALSE, sep = "\t")
    dge_genes = dge_features$V2
    dge_bcs = fread(dge_file("barcodes.tsv"), header = FALSE)$V1

    #- keep 'Gene Expression' features only (e.g. 10x multiome GEX+ATAC Peaks #30)...
    if(ncol(dge_features) >= 3){
        nb_other = sum(dge_features$V3 != "Gene Expression")
        if(nb_other > 0) message("Discarding ", nb_other, " non 'Gene Expression' features (", paste(unique(dge_features$V3[dge_features$V3 != "Gene Expression"]), collapse = ", "), ")")
        dge_genes[dge_features$V3 != "Gene Expression"] = NA_character_
    }
    dge = data.table(
        gene_idx = dge_counts@i + 1L,
        bc_idx = dge_counts@j + 1L,
        counts = dge_counts@x
    )
} else {
    #- Drop-seq tools DGE (GENE + one column per cell barcode), converted column by column...
    dge_counts = fread(args$dge)
    dge_genes = dge_counts[[1]]
    dge_bcs = colnames(dge_counts)[-1]
    dge = rbindlist(lapply(seq_along(dge_bcs), function(j){
        bc_counts = dge_counts[[j + 1]]
        nz = which(bc_counts != 0)
        return(data.table(gene_idx = nz, bc_idx = j, counts = as.numeric(bc_counts[nz])))
    }))
}
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
