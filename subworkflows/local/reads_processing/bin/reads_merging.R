#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(GenomicRanges, quietly = TRUE)
library(dplyr)
options(width = 600)

parser <- ArgumentParser(description = "Merge reads with genomic annotations")
parser$add_argument("--reads", required = TRUE, help = "Path to reads BED file")
parser$add_argument("--annotations", required = TRUE, help = "Path to annotations BED file")
parser$add_argument("--output", required = TRUE, help = "Path to output file")
parser$add_argument("--threads", type = "integer", default = 1, help = "Number of threads")

args <- parser$parse_args()

cat("Reading files...\n")
reads <- fread(args$reads, col.names = c("seqnames", "start", "end", "strand", "read_id", "tags"))
annotations <- fread(args$annotations)[,-15]

cat("Reads:", nrow(reads), "rows\n")
cat("Annotations:", nrow(annotations), "rows\n")

cat("Converting to GRanges objects...\n")

# Convert reads to GRanges
gr_reads <- GRanges(
  seqnames = reads$seqnames,
  ranges = IRanges(start = reads$start, end = reads$end),
  strand = reads$strand,
  read_id = reads$read_id,
  tags = reads$tags
)

# Convert annotations to GRanges
gr_anno <- GRanges(
  seqnames = annotations$seqnames,
  ranges = IRanges(start = annotations$start, end = annotations$end),
  strand = annotations$strand
)

# Add all annotation columns as metadata
anno_mcols <- annotations[, -c("seqnames", "start", "end"), with = FALSE]
mcols(gr_anno) <- anno_mcols

cat("Performing overlap with strand concordance...\n")

# Find overlaps: ignore.strand=FALSE ensures strand must match
overlaps <- findOverlaps(gr_reads, gr_anno, type = "any", ignore.strand = FALSE)

cat("Found", length(overlaps), "overlaps\n")

# Extract overlapping reads and annotations
merged_reads <- gr_reads[queryHits(overlaps)]
merged_anno <- gr_anno[subjectHits(overlaps)]

cat("Converting results to data.table...\n")

# Get the length of results
n_overlaps <- length(overlaps)

# Extract all components
merged <- data.table(
  seqnames = as.character(seqnames(merged_reads)),
  read_start = start(merged_reads),
  read_end = end(merged_reads),
  strand = as.character(strand(merged_reads)),
  read_id = merged_reads$read_id,
  tags = merged_reads$tags,
  anno_start = start(merged_anno),
  anno_end = end(merged_anno)
)

# Add metadata columns individually to avoid issues
anno_mcols_df <- as.data.frame(mcols(merged_anno))
for (col in names(anno_mcols_df)) {
  merged[[col]] <- anno_mcols_df[[col]]
}
print((merged))

#discard suprious mapping
merged %>%
    group_by(tags) %>%
    summarise(
        filt_lastexon = ifelse(
            (exon_number == 1 & strand == "+" & read_start <= anno_start) |
            (exon_number == max(exon_number) & strand == "-" & read_end >= anno_end),
            TRUE, FALSE),
        filt_internalexon = ifelse(
            (exon_number != 1 & strand == "+" & read_start > anno_start),
            TRUE, FALSE)
    ) %>% print()


#writing
fwrite(merged, args$output, sep = "\t", quote = FALSE)