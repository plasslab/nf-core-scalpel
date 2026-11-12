#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))

args = commandArgs(TRUE)
SCRIPT = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
#- ip_filtering.R writes <chr>.ipdb in the working directory...
setwd(OUT)

#- isoforms: X (+, exon 1000-2000, 3' end 2000) ; Y (-, exon 5000-6000, 3' end 5000)
annots = data.table(seqnames = "chr1", start = c(1000, 5000), end = c(2000, 6000), strand = c("+", "-"),
                    startR = 0, endR = 1000, gene_name = c("GX", "GY"), transcript_name = c("X", "Y"), exon_number = 1, TPM_perc = 1)
fwrite(annots, file.path(OUT, "annots.tsv"), sep = "\t")
#- IP sites: X at 1500 (490 nt from end), 1200 (790 nt), 1980 (10 nt: ignored, < threshold) ; Y at 5500 (500 nt)
fwrite(data.table(c("chr1", "chr1", "chr1", "chr1"), c(1500, 1200, 1980, 5500), c(1510, 1210, 1990, 5510), ".", ".", c("+", "+", "+", "-")),
       file.path(OUT, "ip.txt"), sep = "\t", col.names = FALSE)
#- fragments (one read per row)
reads = data.table(
  frag = c("x_near", "x_up", "x_span", "x_span", "y_near", "y_up"),
  read_start = c(1900, 1300, 1400, 1900, 5100, 5700), read_end = c(1997, 1397, 1497, 1997, 5197, 5797),
  transcript_name = c("X", "X", "X", "X", "Y", "Y"), strand = c("+", "+", "+", "+", "-", "-"))
reads[, `:=`(seqnames = "chr1", tags = paste0("CB:Z:", frag, "::UB:Z:A"), gene_name = paste0("G", transcript_name), read_id = paste0("r", .I))]
write_fst(reads, file.path(OUT, "reads.fst"))

run = function(threshold){
  out = file.path(OUT, paste0("out_", threshold, ".fst"))
  status = system2(file.path(R.home("bin"), "Rscript"), c(SCRIPT, "--reads", file.path(OUT, "reads.fst"), "--ip_reference", file.path(OUT, "ip.txt"),
    "--annotations", file.path(OUT, "annots.tsv"), "--output_path", out, "--threshold_dist", threshold, "--threads", "1"),
    stdout = file.path(OUT, "log.txt"), stderr = file.path(OUT, "log.txt"))
  if(status != 0) { cat(readLines(file.path(OUT, "log.txt")), sep = "\n"); quit(status = 1) }
  unique(sub("CB:Z:(.*)::.*", "\\1", read_fst(out)$tags))
}
kept60 = run(60)
kept600 = run(600)
checks = data.table(
  case = c("+ fragment near the 3' end", "+ fragment upstream of the closest IP site (490 nt)", "+ fragment spanning the IP site",
           "- fragment near the 3' end", "- fragment upstream of the IP site", "--threshold_dist 600: closest site now 790 nt, + fragment at 1300 kept"),
  should_keep = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
  kept = c("x_near" %in% kept60, "x_up" %in% kept60, "x_span" %in% kept60, "y_near" %in% kept60, "y_up" %in% kept60, "x_up" %in% kept600))
checks[, result := ifelse(kept == should_keep, "PASS", "FAIL")]
print(checks, row.names = FALSE)
cat(sum(checks$result == "PASS"), "/", nrow(checks), "checks passed\n")
fwrite(checks[, .(case, result)], file.path(OUT, "checks.tsv"), sep = "\t")
