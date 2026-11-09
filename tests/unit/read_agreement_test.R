#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))

args = commandArgs(TRUE)
SCRIPT = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


# Annotation (+ strand): one gene per case, isoforms X (agrees) & Y (one read disagrees)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- exons listed 5'->3' (genomic order on + strand), 1-based inclusive
models = list(
  b_boundary = list(X = list(c(1000, 2000)),                             Y = list(c(1000, 1400), c(1600, 2000))),
  d_junction = list(X = list(c(40000, 40300), c(40500, 41000)),          Y = list(c(40000, 40350), c(40500, 41000))),
  d_exonskip = list(X = list(c(30000, 30200), c(30400, 30600), c(30800, 31000)), Y = list(c(30000, 30200), c(30800, 31000))),
  f6_distance = list(X = list(c(20000, 21000)),                          Y = list(c(20000, 21500)))
)
annots = rbindlist(lapply(names(models), function(g){
  rbindlist(lapply(names(models[[g]]), function(iso){
    ex = rev(models[[g]][[iso]])                        #- 3'->5': exon_number 1 = terminal exon
    len = sapply(ex, function(e) e[2] - e[1] + 1)
    startR = c(0, cumsum(len)[-length(len)])
    data.table(seqnames = "chr1", start = sapply(ex, `[`, 1), end = sapply(ex, `[`, 2), strand = "+",
               startR = startR, endR = startR + len - 1, gene_name = g, transcript_name = paste0(g, "_", iso),
               exon_number = seq_along(ex), TPM_perc = 0.5)
  }))
}))
fwrite(annots, file.path(OUT, "annots.tsv"), sep = "\t")


# Reads: fragment = cell barcode + UMI ; segments of a spliced read share their read id
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- (fragment, read id, start, end) 1-based inclusive
reads = rbindlist(list(
  #- b: read2 crosses Y's exon end into its intron (filter b / F1)
  data.table(frag = "b1", rid = "b1r1", start = 1700, end = 1797),
  data.table(frag = "b1", rid = "b1r2", start = 1350, end = 1447),
  #- control b: both reads inside both isoforms' terminal exon
  data.table(frag = "b2", rid = "b2r1", start = 1700, end = 1797),
  data.table(frag = "b2", rid = "b2r2", start = 1800, end = 1897),
  #- d junction: read2 spliced at X's donor (40300), inside Y's longer exon (filter d / F3)
  data.table(frag = "j1", rid = "j1r1", start = 40900, end = 40997),
  data.table(frag = "j1", rid = c("j1r2", "j1r2"), start = c(40250, 40500), end = c(40300, 40547)),
  #- d exon skip: read2 spliced exon2->exon1 of X, skipping nothing in X but Y lacks exon 30400-30600 (filter d / F2)
  data.table(frag = "e1", rid = "e1r1", start = 30900, end = 30997),
  data.table(frag = "e1", rid = c("e1r2", "e1r2"), start = c(30550, 30800), end = c(30600, 30847)),
  #- F6: read2 is 120 nt from X's end but 620 nt from Y's end (> 600)
  data.table(frag = "f1", rid = "f1r1", start = 20893, end = 20990),
  data.table(frag = "f1", rid = "f1r2", start = 20783, end = 20880),
  #- control F6: both reads within 600 nt of both ends
  data.table(frag = "f2", rid = "f2r1", start = 20903, end = 21000),
  data.table(frag = "f2", rid = "f2r2", start = 20943, end = 20990)
))
bed = data.table(V1 = "chr1", V2 = reads$start - 1, V3 = reads$end, V4 = reads$rid, V5 = 255, V6 = "+")
for(i in 7:13) bed[[paste0("V", i)]] = "."
bed$V14 = paste0("CB:Z:", reads$frag)
bed$V15 = "UB:Z:AAAA"
fwrite(bed, file.path(OUT, "reads.bed"), sep = "\t", col.names = FALSE)

status = system2(file.path(R.home("bin"), "Rscript"), c(SCRIPT, "--reads", file.path(OUT, "reads.bed"), "--annotations", file.path(OUT, "annots.tsv"),
  "--output", file.path(OUT, "out.fst"), "--distance_threshold", "600"), stdout = file.path(OUT, "log.txt"), stderr = file.path(OUT, "log.txt"))
if(status != 0) { cat(readLines(file.path(OUT, "log.txt")), sep = "\n"); quit(status = 1) }


# Expected fragment x isoform pairs kept
# ++++++++++++++++++++++++++++++++++++++
kept = unique(read_fst(file.path(OUT, "out.fst"), as.data.table = TRUE)[, .(frag = sub("CB:Z:(.*)::.*", "\\1", tags), transcript_name)])
expected = data.table(
  frag = c("b1", "b1", "b2", "b2", "j1", "j1", "e1", "e1", "f1", "f1", "f2", "f2"),
  isoform = c("X", "Y", "X", "Y", "X", "Y", "X", "Y", "X", "Y", "X", "Y"),
  case = c("b: read crosses Y exon boundary", "b: read crosses Y exon boundary", "control b", "control b",
           "d: read spliced at another donor than Y", "d: read spliced at another donor than Y",
           "d: spliced read on exons absent from Y", "d: spliced read on exons absent from Y",
           "F6: one read > 600 nt from Y end", "F6: one read > 600 nt from Y end", "control F6", "control F6"),
  should_keep = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
)
gene_of = c(b = "b_boundary", j = "d_junction", e = "d_exonskip", f = "f6_distance")
expected[, transcript_name := paste0(gene_of[substr(frag, 1, 1)], "_", isoform)]
expected[, kept := paste(frag, transcript_name) %in% paste(kept$frag, kept$transcript_name)]
expected[, result := ifelse(kept == should_keep, "PASS", "FAIL")]
print(expected[, .(case, frag, isoform, should_keep, kept, result)], row.names = FALSE)
cat(sum(expected$result == "PASS"), "/", nrow(expected), "checks passed\n")
fwrite(expected[, .(case, result)], file.path(OUT, "checks.tsv"), sep = "\t")
