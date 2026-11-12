#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))


# Argument Parser
# +++++++++++++++
parser = ArgumentParser(description='Internal priming reference generation (genomic A-rich regions)')
parser$add_argument('--sequence', type="character", help='path of chromosome sequence (single line, no header)')
parser$add_argument('--seqname', type="character", help='chromosome name')
parser$add_argument('--output', type="character", help='path of output internal priming sites file')
parser$add_argument('--window', type="character", default="10", help='window size (nt)')
parser$add_argument('--min_a', type="character", default="7", help='minimum number of A in a window')
parser$add_argument('--min_run', type="character", default="6", help='minimum stretch of consecutive A starting a window')
args = parser$parse_args()

WINDOW = as.integer(args$window)
MIN_A = as.integer(args$min_a)
MIN_RUN = as.integer(args$min_run)
CHUNK = 10000000L


# Functions
# +++++++++
#- windows (1-based start) with >= MIN_A bases, or a stretch of MIN_RUN bases at their 3' end (+: start, -: end)...
window_hits = function(chr_seq, base, strand){
  nb_windows = length(chr_seq) - WINDOW + 1
  if(nb_windows < 1) return(integer(0))
  hits = lapply(seq(1, nb_windows, by = CHUNK), function(from){
    to = min(from + CHUNK - 1, nb_windows)
    cs = c(0L, cumsum(chr_seq[from:(to + WINDOW - 1)] == base))
    j = seq_len(to - from + 1)
    nb_base = cs[j + WINDOW] - cs[j]
    if(strand == "+"){
      run = cs[j + MIN_RUN] - cs[j]
    } else {
      run = cs[j + WINDOW] - cs[j + WINDOW - MIN_RUN]
    }
    return(as.integer(from) - 1L + which(nb_base >= MIN_A | run == MIN_RUN))
  })
  return(unlist(hits))
}

#- merge overlapping windows into sites & trim their 5' end to the first A (+) / last T (-)...
merge_windows = function(chr_seq, hits, base, strand){
  if(length(hits) == 0) return(data.table())
  new_site = c(TRUE, diff(hits) >= WINDOW)
  sites = data.table(start = hits[new_site], end = hits[c(new_site[-1], TRUE)] + WINDOW - 1L)
  for(k in seq_len(WINDOW - MIN_RUN + 1)){
    if(strand == "+"){
      sites[chr_seq[start] != base, start := start + 1L]
    } else {
      sites[chr_seq[end] != base, end := end - 1L]
    }
  }
  return(sites[, .(seqnames = args$seqname, start = start - 1L, end, V4 = ".", V5 = ".", strand = strand)])
}


# Internal priming sites
# ++++++++++++++++++++++
message("Reading ", args$seqname, " sequence...")
chr_seq = readBin(args$sequence, what = "raw", n = file.size(args$sequence))
#- upper case (soft-masked genomes)...
chr_seq = chr_seq & as.raw(0xDF)

message("Scanning A-rich windows...")
ip_sites = rbindlist(list(
  merge_windows(chr_seq, window_hits(chr_seq, charToRaw("A"), "+"), charToRaw("A"), "+"),
  merge_windows(chr_seq, window_hits(chr_seq, charToRaw("T"), "-"), charToRaw("T"), "-")
))
message(args$seqname, ": ", nrow(ip_sites), " internal priming sites")

fwrite(ip_sites, file = args$output, sep = "\t", col.names = FALSE)
