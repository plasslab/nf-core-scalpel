#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))

args = commandArgs(TRUE)
SCRIPT = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


# Synthetic chromosome: one case per 40 nt block (C background, no A/T)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cases = c(
  count_7a    = "AACAAGAACA",                     #- 10 nt with 7 A: site
  count_6a    = "AACAAGACCA",                     #- 10 nt with 6 A, no stretch: no site
  run_6a      = "AAAAAAGGCC",                     #- stretch of 6 A (6 A only): site AAAAAAGGCC
  trim_5p     = "GAAACAAGAA",                     #- windows at +0 & +1 (7 A) merged, trimmed to the first A
  bookended   = "AAAACAACAGCAACAACAAAAACAG",      #- bookended windows: 2 sites
  minus_7t    = "TGTTCTTGTT",                     #- 7 T in 10 nt: - strand site
  minus_trim  = "TTGTTCTTTG",                     #- windows at -1 & +0 (7 T) merged, trimmed to the last T (- strand 5' end)
  lowercase   = "aacaagaaca"                      #- soft-masked: site
)
block = 40
chr_seq = paste0(sapply(cases, function(x) paste0(x, strrep("C", block - nchar(x)))), collapse = "")
writeChar(chr_seq, file.path(OUT, "chrT.seq"), eos = NULL)
starts = (seq_along(cases) - 1) * block

status = system2(file.path(R.home("bin"), "Rscript"), c(SCRIPT, "--sequence", file.path(OUT, "chrT.seq"), "--seqname", "chrT",
  "--output", file.path(OUT, "chrT.ip.txt")), stdout = file.path(OUT, "log.txt"), stderr = file.path(OUT, "log.txt"))
if(status != 0) { cat(readLines(file.path(OUT, "log.txt")), sep = "\n"); quit(status = 1) }
sites = fread(file.path(OUT, "chrT.ip.txt"), header = FALSE, col.names = c("seqnames", "start", "end", "V4", "V5", "strand"))
site = function(from, to, strand) paste(from, to, strand) %in% paste(sites$start, sites$end, sites$strand)
in_block = function(i) sites[start >= starts[i] & end <= starts[i] + block]


# Expected sites (0-based start, end)
# +++++++++++++++++++++++++++++++++++
checks = data.table(
  case = c("10 nt window with 7 A -> site", "10 nt window with 6 A only -> no site", "stretch of 6 A -> site",
           "site trimmed to its first A", "bookended windows -> 2 sites", "10 nt window with 7 T -> - strand site",
           "- strand site trimmed to its last T", "soft-masked (lower case) sequence -> site",
           "integer coordinates (no scientific notation)"),
  pass = c(
    site(starts[1], starts[1] + 10, "+"),
    nrow(in_block(2)) == 0,
    site(starts[3], starts[3] + 10, "+"),
    site(starts[4] + 1, starts[4] + 11, "+"),
    site(starts[5], starts[5] + 10, "+") & site(starts[5] + 11, starts[5] + 25, "+"),
    site(starts[6], starts[6] + 10, "-"),
    site(starts[7] - 1, starts[7] + 9, "-"),
    site(starts[8], starts[8] + 10, "+"),
    !any(grepl("e", readLines(file.path(OUT, "chrT.ip.txt"))))
  ))
checks[, result := ifelse(pass, "PASS", "FAIL")]
print(checks[, .(case, result)], row.names = FALSE)
cat(sum(checks$result == "PASS"), "/", nrow(checks), "checks passed\n")
fwrite(checks[, .(case, result)], file.path(OUT, "checks.tsv"), sep = "\t")
