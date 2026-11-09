#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse,   quietly = TRUE)
  library(data.table, quietly = TRUE)
  library(dplyr,     quietly = TRUE)
  library(fst,        quietly = TRUE)
})

parser <- ArgumentParser(description = "Merge reads with genomic annotations")
parser$add_argument("--reads", required = TRUE, help = "Path to reads BED file")
parser$add_argument("--annotations", required = TRUE, help = "Path to annotations BED file")
parser$add_argument("--output", required = TRUE, help = "Path to output file")
parser$add_argument("--distance_threshold", required = TRUE, type = "integer",
  help = "Transcriptomic distance threshold for discarding reads")
parser$add_argument("--threads", type = "integer", default = 1L, help = "Number of threads")
args <- parser$parse_args()

setDTthreads(args$threads)

print(args$reads)

# ---------------------------------------------------------------------------
# 1. Read input files
# ---------------------------------------------------------------------------
message("step1...")
reads = fread(
  file = args$reads,
  select = c(1, 2, 3, 6, 4, 14, 15),
  col.names = c("seqnames", "start", "end", "strand", "read_id", "bc", "umi"),
  nThread = args$threads,
  fill = TRUE
)
reads = reads %>%
  dplyr::mutate(tags = paste(bc, umi, sep = "::")) %>%
  dplyr::select(-bc, -umi)
reads[, seg_id := .I]
reads[, read_id := sub("/.*", "", read_id)]
reads[, c("nsplices", "start") := .(uniqueN(start), start + 1L), by = read_id]

annotations = fread(
  args$annotations,
  select = c("seqnames", "start", "end", "strand",
             "startR", "endR", "gene_name", "transcript_name",
             "exon_number", "TPM_perc"),
  nThread = args$threads
)

# ---------------------------------------------------------------------------
# 2. Strand-specific overlap using foverlaps (replaces GenomicRanges)
#    Key annotations on (seqnames, strand, start, end).
#    foverlaps matches seqnames+strand exactly, then does interval overlap.
# ---------------------------------------------------------------------------
message("step2...")

# Rename annotation coords to avoid collision (foverlaps adds both)
setnames(annotations, c("start", "end"), c("anno_start", "anno_end"))

# Key annotations — last two key columns must be the interval
setkey(annotations, seqnames, strand, anno_start, anno_end)

# Reads need matching key columns; foverlaps uses start/end from x by default
# but we need to tell it which columns in x are the interval
merged <- foverlaps(
  reads, annotations,
  by.x = c("seqnames", "strand", "start", "end"),
  by.y = c("seqnames", "strand", "anno_start", "anno_end"),
  type = "any",
  nomatch = NULL
)

# ---- Filter a — discard fragments associated to intergenic/intronic reads (#37) ----
#  (any read segment without strand-specific exon overlap discards the whole fragment)
unmapped_tags <- reads[!merged, on = "seg_id", unique(tags)]
merged <- merged[!(tags %in% unmapped_tags)]
rm(reads, annotations, unmapped_tags)

# Rename read coords to match v2 output naming
setnames(merged, c("start", "end"), c("read_start", "read_end"))

# ---------------------------------------------------------------------------
# 3. Integer codes for grouping
# ---------------------------------------------------------------------------
message("step3...")
merged[, tags_i       := .GRP, by = tags]
merged[, read_id_i    := .GRP, by = read_id]
merged[, transcript_i := .GRP, by = transcript_name]
merged[, fragTr_id    := .GRP, by = .(tags_i, transcript_i)]
merged[, readTr_id    := .GRP, by = .(read_id_i, transcript_i)]

# ---------------------------------------------------------------------------
# 4. Filter chain
#    Same order as original mapping_filtering.R: a (above) -> b (F1) -> c (F4)
#    -> d (F3, F2, F5). F4 must see reads after F1 only (#37).
#    F6 then F7 sequential (F7 depends on F6)
# ---------------------------------------------------------------------------
message("step4...")
# ---- Filter 1 — reads overlapping transcript coordinates (#b) ----
merged[, f_bad := (
  (exon_number == 1L & strand == "+" & read_start < anno_start) |
  (exon_number == 1L & strand == "-" & read_end   > anno_end)   |
  (exon_number != 1L & (read_start < anno_start | read_end > anno_end))
)]
bad <- data.table(fragTr_id = merged[f_bad == TRUE, unique(fragTr_id)])
merged[, f_bad := NULL]
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---- Filter 4 — fragment completeness (#c, #37) ----
#  keep fragment-transcript pairs only if the transcript intersects ALL the
#  read ids of the fragment (nb_reads_counted == nb_reads_tot)
merged[, nb_reads_tot := uniqueN(read_id_i), by = tags_i]
tmp <- merged[, .(nb_reads_counted = uniqueN(read_id_i), nb_reads_tot = nb_reads_tot[1L]), by = fragTr_id]
bad <- data.table(fragTr_id = tmp[nb_reads_counted != nb_reads_tot, unique(fragTr_id)])
merged[, nb_reads_tot := NULL]
rm(tmp)
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---- Filter 3 — spliced reads coherency (#d) ----
bad <- data.table(fragTr_id = merged[
  nsplices > 1L & read_start != anno_start & read_end != anno_end,
  unique(fragTr_id)
])
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---- Filter 2 — exon consecutivity ----
bad_rids <- merged[
  nsplices > 1L,
  .(exon_rng = max(exon_number) - min(exon_number) + 1L, ns = nsplices[1L]),
  by = readTr_id
][exon_rng != ns, .(readTr_id)]
bad <- data.table(fragTr_id = unique(merged[bad_rids, fragTr_id, on = "readTr_id", nomatch = NULL]))
rm(bad_rids)
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---- Filter 5 — spliced/unspliced concordance ----
spliced_tags  <- data.table(tags_i = merged[nsplices > 1L, unique(tags_i)])
spliced_frags <- unique(merged[nsplices > 1L, .(tags_i, fragTr_id)])
candidates <- merged[nsplices == 1L][spliced_tags, on = "tags_i", nomatch = NULL]
bad <- data.table(fragTr_id = candidates[!spliced_frags, on = .(tags_i, fragTr_id), unique(fragTr_id)])
rm(spliced_tags, spliced_frags, candidates)
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---------------------------------------------------------------------------
# 5. Relative transcriptomic coordinates
# ---------------------------------------------------------------------------
message("step5...")
merged[strand == "+", `:=`(
  start.rdR = endR - (read_end   - anno_start),
  end.rdR   = endR - (read_start - anno_start)
)]
merged[strand == "-", `:=`(
  start.rdR = endR - (anno_end - read_start),
  end.rdR   = endR - (anno_end - read_end)
)]

# ---- Filter 6 — distance threshold ----
message("step6...")
bad <- data.table(fragTr_id = merged[start.rdR > args$distance_threshold, unique(fragTr_id)])
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

# ---- Filter 7 — reads must reach last exon ----
message("step7...")
bad_tr <- merged[, .(min_exon = min(exon_number)), by = transcript_i
                 ][min_exon != 1L, .(transcript_i)]
bad <- data.table(fragTr_id = unique(merged[bad_tr, fragTr_id, on = "transcript_i", nomatch = NULL]))
rm(bad_tr)
merged <- merged[!bad, on = "fragTr_id"]
rm(bad)

invisible(gc())

# ---------------------------------------------------------------------------
# 6. Finalise and write output
# ---------------------------------------------------------------------------
message("step8...")
merged[, fg_start.rdR := min(start.rdR), by = fragTr_id]
merged[, fg_end.rdR   := max(end.rdR),   by = fragTr_id]

merged[, c("nsplices", "fragTr_id", "readTr_id", "seg_id",
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

fst::write_fst(merged, path = args$output)
