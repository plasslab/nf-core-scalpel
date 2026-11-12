#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))

args = commandArgs(TRUE)
SCRIPT = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
#- isoform_selection.R writes stats.txt in the working directory...
setwd(OUT)
RSCRIPT = file.path(R.home("bin"), "Rscript")


# Annotation (+ strand), exons listed 5' -> 3' (genomic order)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
models = list(
  #- same structure, 3' ends 30 nt apart -> collapsed (B, highest TPM, extended to A's 3' end)...
  G1 = list(A = list(c(500, 700), c(1000, 2000)), B = list(c(500, 700), c(1000, 1970))),
  #- same terminal exon, upstream exon of the same length at another position (inside the 600 nt window) -> kept...
  G2 = list(C = list(c(5000, 5100), c(6000, 6300)), D = list(c(5400, 5500), c(6000, 6300))),
  #- differences only beyond 600 nt from the 3' end -> collapsed...
  G3 = list(E = list(c(9000, 9100), c(9500, 9600), c(10000, 10700)), F = list(c(8000, 8100), c(9500, 9600), c(10000, 10700))),
  #- same structure, 3' ends 80 nt apart (> --distance_3end 50) -> kept...
  G4 = list(G = list(c(20000, 21000)), H = list(c(20000, 21080)))
)
tpm = c(A = 0.4, B = 0.6, C = 0.5, D = 0.5, E = 0.7, F = 0.3, G = 0.5, H = 0.5)
tx = rbindlist(lapply(names(models), function(g){
  rbindlist(lapply(names(models[[g]]), function(iso){
    ex = rev(models[[g]][[iso]])                         #- 3' -> 5': exon_number 1 = terminal exon
    w = sapply(ex, function(e) e[2] - e[1])
    startR = endR = numeric(length(ex))
    for(i in seq_along(ex)){
      startR[i] = if(i == 1) 0 else endR[i - 1] + 1
      endR[i] = if(i == 1) w[i] - 1 else startR[i] + w[i]
    }
    data.table(seqnames = "chrT", start = sapply(ex, `[`, 1), end = sapply(ex, `[`, 2), width = w, strand = "+",
               startR = startR, endR = endR, gene_id = g, gene_name = g, transcript_id = iso, transcript_name = iso,
               exon_number = seq_along(ex), TPM_perc = tpm[[iso]])
  }))
}))
fwrite(tx, file.path(OUT, "chrT.tsv"), sep = "\t")

status = system2(RSCRIPT, c(SCRIPT, "--input", file.path(OUT, "chrT.tsv"), "--output", file.path(OUT, "chrT.annots.gz"),
  "--distance_3end", "50", "--distance_profile", "600"), stdout = file.path(OUT, "log.txt"), stderr = file.path(OUT, "log.txt"))
if(status != 0) { cat(readLines(file.path(OUT, "log.txt")), sep = "\n"); quit(status = 1) }
annots = fread(file.path(OUT, "chrT.annots.gz"))
kept = unique(annots$transcript_name)
group_of = function(iso) unique(annots[transcript_name == iso]$transcript_collapseds)


# Checks
# ++++++
checks = data.table(
  case = c("same structure, 3' ends 30 nt apart -> collapsed, highest TPM kept", "representative extended to the furthest 3' end of its group",
           "relative coordinates shifted with the extension", "upstream exon of the same length at another position -> not collapsed",
           "differences beyond 600 nt only -> collapsed", "3' ends 80 nt apart (> distance_3end) -> not collapsed"),
  pass = c(
    identical(sort(intersect(kept, c("A", "B"))), "B") && identical(group_of("B"), "A,B"),
    annots[transcript_name == "B" & startR == 0]$end == 2000,
    annots[transcript_name == "B" & startR == 0]$endR == 999 && annots[transcript_name == "B" & startR > 0]$startR == 1000,
    all(c("C", "D") %in% kept),
    identical(sort(intersect(kept, c("E", "F"))), "E") && identical(group_of("E"), "E,F"),
    all(c("G", "H") %in% kept)
  ))
checks[, result := ifelse(pass, "PASS", "FAIL")]
print(checks[, .(case, result)], row.names = FALSE)
cat(sum(checks$result == "PASS"), "/", nrow(checks), "checks passed\n")
fwrite(checks[, .(case, result)], file.path(OUT, "checks.tsv"), sep = "\t")
