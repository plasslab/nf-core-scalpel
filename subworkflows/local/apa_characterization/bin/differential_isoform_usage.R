#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(Seurat, quietly = TRUE))
source("differential_isoform_usage_functions.R") %>% suppressPackageStartupMessages()


#ArgumentParser
#==============
parser = ArgumentParser(description='Differential isoform usage between clusters')
parser$add_argument('--idge', type="character", help='path of directory containing the <sample>_iDGE directories')
parser$add_argument('--clusters', type="character", default="NULL", help='path of clusters file (orig.ident, Barcodes, Clusters) or NULL')
parser$add_argument('--output', type="character", help='path of output DIU table')
parser$add_argument('--output_rds', type="character", help='path of output seurat file')
parser$add_argument('--annotations', type="character", default="NULL", help='path of directory containing the isoform annotation (*.annots.gz) or NULL')
parser$add_argument('--test', type="character", default="chisq", help='DIU test: chisq (default, original) or permutation')
parser$add_argument('--min_delta', type="double", default=0.1, help='minimal change of isoform usage between groups (default: 0.1)')
parser$add_argument('--permutations', type="integer", default=10000, help='maximal number of cell permutations (test permutation, default: 10000)')
parser$add_argument('--cpus', type="integer", default=1, help='number of cores')
args = parser$parse_args()


THR_ABUND = 0.05
THR_PVAL = 0.05
THR_VAR = args$min_delta


#get iDGE directories
all_dirs = list.files(path=args$idge, pattern = "_iDGE$", full.names = T)
if(length(all_dirs) == 0) stop("No <sample>_iDGE directory found in ", args$idge)

#Opening (1)
#*******
all.objs = lapply(all_dirs, function(x){
  sample_name = sub("_iDGE$", "", basename(x))
  message("Loading ", sample_name, "...")
  counts = Seurat::Read10X(data.dir = x, gene.column = 1)
  curr.obj = Seurat::CreateSeuratObject(counts, project = sample_name, min.cells = 0, min.features = 0)
  #- explicit sample & barcode annotation (independent of cell renaming at merge)...
  curr.obj$orig.ident = sample_name
  curr.obj$Barcodes = colnames(curr.obj)
  return(curr.obj)
})
names(all.objs) = sub("_iDGE$", "", basename(all_dirs))

#merging (cell names prefixed by sample to avoid barcode collisions)
if(length(all.objs)>1){
	seurat.obj = merge(all.objs[[1]], all.objs[2:length(all.objs)], add.cell.ids = names(all.objs))
}else{
	seurat.obj = all.objs[[1]]
}

#annotate cell barcodes with clusters if needed (default samples) (2)
#****************************************************************
if(args$clusters == "NULL"){
	seurat.obj$Clusters = seurat.obj$orig.ident
}else{
	#a. read clusters table (3 columns: orig.ident = sample name, Barcodes = cell barcode, Clusters)
	cluster_df = fread(args$clusters, header = "auto")
	if(ncol(cluster_df) != 3){
		stop("Clusters file must have exactly 3 columns (sample, barcode, cluster), found ", ncol(cluster_df), ": ", paste(colnames(cluster_df), collapse = ", "))
	}
	colnames(cluster_df) = c("orig.ident", "Barcodes", "Clusters")
	cluster_df = cluster_df %>%
		dplyr::mutate(across(everything(), as.character)) %>%
		dplyr::distinct(orig.ident, Barcodes, .keep_all = T)

	#b merge cluster into the seurat object
	meta_tab = seurat.obj@meta.data %>%
		dplyr::select(orig.ident, Barcodes) %>%
		dplyr::left_join(cluster_df, by = c("orig.ident", "Barcodes"))
	nb_matched = sum(!is.na(meta_tab$Clusters))
	message("Clusters: ", nb_matched, " / ", nrow(meta_tab), " cells annotated")
	if(nb_matched == 0){
		stop("No cell of the clusters file matches the iDGE cells. Expected sample names: ",
			paste(unique(seurat.obj$orig.ident), collapse = ", "), " ; barcodes like: ",
			paste(head(seurat.obj$Barcodes, 3), collapse = ", "), " (got e.g. ",
			paste(head(cluster_df$orig.ident, 1), head(cluster_df$Barcodes, 1)), ")")
	}
	seurat.obj$Clusters = meta_tab$Clusters
	#- cells without cluster are not used in the differential analysis...
	if(nb_matched < nrow(meta_tab)){
		message("Discarding ", nrow(meta_tab) - nb_matched, " cells without cluster annotation")
		seurat.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(seurat.obj$Clusters)])
	}
}

print(seurat.obj)

#Differential analysis (3)
#*********************
if(n_distinct(seurat.obj$Clusters) == 1){
	my_isoforms = data.table(output="Provide more condition with levels > 1 for differential isoform usage analysis")
}else{
	#- isoform 3' ends (direction of the 3' end shifts)...
	annot_files = if(args$annotations == "NULL") character(0) else list.files(args$annotations, pattern = "annots.gz$", full.names = T)
	ends = Read_isoformEnds(annot_files)
	#- cells permuted within samples when clusters are compared...
	strata = if(args$clusters == "NULL") NULL else seurat.obj$orig.ident
	my_isoforms = FindIsoforms(
		counts = SeuratObject::LayerData(JoinLayers(seurat.obj), layer = "counts"),
		groups = seurat.obj$Clusters,
		strata = strata,
		ends = ends,
		test = args$test,
		threshold_pval = THR_PVAL,
		threshold_abund = THR_ABUND,
		threshold_var = THR_VAR,
		n_perm = args$permutations,
		cores = args$cpus)
}

#Writing differential analysis table (4)
#***********************************
fwrite(my_isoforms, file=args$output, sep="\t", nThread=1)
saveRDS(seurat.obj, file=args$output_rds)
