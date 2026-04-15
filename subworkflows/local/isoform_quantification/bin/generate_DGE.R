
suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(Seurat, quietly = TRUE))
suppressPackageStartupMessages(library(reshape2, quietly = TRUE))


parser = ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('--dge', type="character", help='path of DGE file')
parser$add_argument('--pred', type="character", help='path of predicted isoforms')
parser$add_argument('--output', type="character", help="path of output file")
args = parser$parse_args()


# File opening

#merge all em files
list_em_files = list.files(args$pred, full.names = F, pattern = "*em_results.fst")

#load cell isoform probabilities
preds = lapply(list_em_files, function(x) {
    tmp = read_fst(x, as.data.table = TRUE)
    return(tmp)
}) %>% rbindlist()

#load DGE count matrix
message("Loading DGE matrix...")
dge = data.frame(Seurat::Read10X(data.dir = args$dge))
colnames(dge) = gsub("\\.", "\\-", colnames(dge))
dge = data.table(dge, keep.rownames = "gene_name")

#melting
message("Melting DGE matrix...")
dge = reshape2::melt(dge)

#renaming
colnames(dge) = c("gene_name","bc","counts")

#subset cells into predictions
dge = dge %>% dplyr::filter(bc %in% preds$bc)

#left_joining by barcodes & gene_name
message("Merging DGE & isoform predictions...")
dge_preds = left_join(preds, dge, by=c("bc", "gene_name")) %>%
    mutate(counts=ifelse(is.na(counts), 0, counts)) %>%
    mutate(tr_counts=as.numeric(rel_abund * counts), gene_transcript = paste0(gene_name, "***", transcript_name)) %>% 
    distinct(gene_transcript, bc, tr_counts) %>%
    tidyr::pivot_wider(names_from=bc, values_from=tr_counts, values_fill=0)

#writing
fwrite(dge_preds, file=args$output, sep="\t")