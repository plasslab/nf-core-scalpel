#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse,      quietly = TRUE)
  library(data.table,    quietly = TRUE)
  library(GenomicRanges, quietly = TRUE)
  library(fst, quietly = TRUE)
})
options(width = 600)

parser <- ArgumentParser(description = "Merge reads with genomic annotations")
parser$add_argument("--reads", required = TRUE, help = "Path to reads BED file")
parser$add_argument("--annotations", required = TRUE, help = "Path to annotations BED file")
parser$add_argument("--output", required = TRUE, help = "Path to output file")
parser$add_argument("--distance_threshold", required = TRUE, type = "integer",
  help = "Transcriptomic distance threshold for discarding reads")
parser$add_argument("--threads", type = "integer", help = "Number of threads")
args <- parser$parse_args()


# 1. Read input files
reads <- fread(
  args$reads,
  col.names = c("seqnames", "start", "end", "strand", "read_id", "tags"),
  nThread = args$threads
)
reads[, read_id := sub("/.*", "", read_id)]
reads[, c("nsplices", "start") := .(uniqueN(start), start + 1L), by = read_id]

annotations <- fread(
  args$annotations,
  select = c("seqnames", "start", "end", "strand",
             "startR", "endR", "gene_name", "transcript_name",
             "exon_number", "TPM_perc"),
  nThread = args$threads
)


# 2. Strand-specific overlap
gr_reads <- GRanges(reads$seqnames, IRanges(reads$start, reads$end),
                    strand = reads$strand)
gr_anno  <- GRanges(annotations$seqnames, IRanges(annotations$start, annotations$end),
                    strand = annotations$strand)

ov <- findOverlaps(gr_reads, gr_anno, type = "any", ignore.strand = FALSE)
rm(gr_reads, gr_anno); invisible(gc())

q <- queryHits(ov)
s <- subjectHits(ov)
rm(ov)


# 3. Build merged table then immediately free source tables
merged <- cbind(
  reads[q,       .(seqnames, read_start = start, read_end = end,
                   strand, read_id, tags, nsplices)],
  annotations[s, .(anno_start = start, anno_end = end, startR, endR,
                   gene_name, transcript_name, exon_number, TPM_perc)]
)
rm(reads, annotations, q, s); invisible(gc())

# 4. Integer codes for grouping (keep original strings for output)
tags_lvls           <- unique(merged$tags)
read_id_lvls        <- unique(merged$read_id)
transcript_lvls     <- unique(merged$transcript_name)

merged[, tags_i        := match(tags,           tags_lvls)]
merged[, read_id_i     := match(read_id,         read_id_lvls)]
merged[, transcript_i  := match(transcript_name, transcript_lvls)]
rm(tags_lvls, read_id_lvls, transcript_lvls); invisible(gc())

# Group IDs using integer codes
merged[, fragTr_id := .GRP, by = .(tags_i, transcript_i)]
merged[, readTr_id := .GRP, by = .(read_id_i, transcript_i)]

# 5. Filter chain

# Filter 1 — boundary containment
bad <- merged[
  (exon_number == 1L & strand == "+" & read_start < anno_start) |
  (exon_number == 1L & strand == "-" & read_end   > anno_end)   |
  (exon_number != 1L & (read_start < anno_start | read_end > anno_end)),
  unique(fragTr_id)]
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# Filter 2 — exon consecutivity
bad_rids <- merged[
  nsplices > 1L,
  .(exon_rng = max(exon_number) - min(exon_number) + 1L, ns = nsplices[1L]),
  by = readTr_id
][exon_rng != ns, readTr_id]
bad <- merged[readTr_id %in% bad_rids, unique(fragTr_id)]
rm(bad_rids)
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# Filter 3 — junction precision
bad <- merged[
  nsplices > 1L & read_start != anno_start & read_end != anno_end,
  unique(fragTr_id)]
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# Filter 4 — fragment completeness (max transcript coverage per tag)
tmp <- merged[, .(n_reads = uniqueN(read_id_i)), by = .(tags_i, transcript_i, fragTr_id)]
tmp[, max_reads := max(n_reads), by = tags_i]
bad <- tmp[n_reads != max_reads, unique(fragTr_id)]
rm(tmp); invisible(gc())
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# Filter 5 — spliced/unspliced concordance
spliced_frags <- unique(merged[nsplices > 1L, .(tags_i, fragTr_id)])
bad <- merged[
  nsplices == 1L &
  tags_i     %in% spliced_frags$tags_i &
  !fragTr_id %in% spliced_frags$fragTr_id,
  unique(fragTr_id)]
rm(spliced_frags); invisible(gc())
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# 6. Relative transcriptomic coordinates
merged[strand == "+", `:=`(
  start.rdR = endR - (read_end   - anno_start),
  end.rdR   = endR - (read_start - anno_start)
)]
merged[strand == "-", `:=`(
  start.rdR = endR - (anno_end - read_start),
  end.rdR   = endR - (anno_end - read_end)
)]

# Filter 6 — distance threshold
bad <- merged[start.rdR > args$distance_threshold, unique(fragTr_id)]
merged <- merged[!fragTr_id %in% bad]; invisible(gc())

# Filter 7 — reads must reach last exon (transcript-level, not per-fragment)
bad_tr <- merged[, .(min_exon = min(exon_number)), by = transcript_i
                 ][min_exon != 1L, unique(transcript_i)]
bad    <- merged[transcript_i %in% bad_tr, unique(fragTr_id)]
rm(bad_tr)
merged <- merged[!fragTr_id %in% bad]; invisible(gc())


# 7. Finalise and write output
merged[, fg_start.rdR := min(start.rdR), by = fragTr_id]
merged[, fg_end.rdR   := max(end.rdR),   by = fragTr_id]

# Remove all helper columns
merged[, c("nsplices", "fragTr_id", "readTr_id",
           "tags_i", "read_id_i", "transcript_i") := NULL]

setorder(merged, transcript_name, fg_start.rdR)
merged <- unique(merged)

setcolorder(merged, c(
  "seqnames", "read_start", "read_end", "strand",
  "anno_start", "anno_end", "startR", "endR",
  "start.rdR", "end.rdR",
  "gene_name", "transcript_name", "TPM_perc", "exon_number",
  "tags", "read_id",
  "fg_start.rdR", "fg_end.rdR"
))

fst::write_fst(
  merged,
  path = args$output
)