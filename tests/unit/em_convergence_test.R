#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))

args = commandArgs(TRUE)
SCRIPT = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RSCRIPT = file.path(R.home("bin"), "Rscript")


# Fragments (one cell): probability of each fragment for each compatible isoform
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- GD (dominated): 10 fragments, all closer to the 3' end of A than of B (B never the best explanation)...
gd = data.table(gene_name = "GD", umi = rep(sprintf("d%02d", 1:10), each = 2), transcript_name = rep(c("GD-A", "GD-B"), 10),
                frag_probs_weighted = rep(c(0.30, 0.25), 10))
#- GS (shared): 6 fragments specific to A, 2 specific to B -> 0.75 / 0.25...
gs = data.table(gene_name = "GS", umi = sprintf("s%02d", 1:8), transcript_name = rep(c("GS-A", "GS-B"), c(6, 2)), frag_probs_weighted = 0.2)
#- GZ: single isoform, fragment beyond the learned distances (probability 0) -> still 1...
gz = data.table(gene_name = "GZ", umi = "z01", transcript_name = "GZ-A", frag_probs_weighted = 0)
cell = rbindlist(list(gd, gs, gz))[, bc := "AAAC-1"]
setcolorder(cell, c("bc", "gene_name", "transcript_name", "umi", "frag_probs_weighted"))
write_fst(cell, file.path(OUT, "cell.fst"))

run = function(extra){
  out = file.path(OUT, paste0("em_", length(extra), ".fst"))
  status = system2(RSCRIPT, c(SCRIPT, "--cell", file.path(OUT, "cell.fst"), "--output_path", out, extra),
    stdout = file.path(OUT, "log.txt"), stderr = file.path(OUT, "log.txt"))
  if(status != 0) { cat(readLines(file.path(OUT, "log.txt")), sep = "\n"); quit(status = 1) }
  em = read_fst(out, as.data.table = TRUE)
  return(setNames(em$rel_abund, em$transcript_name))
}
converged = run(character(0))
original = run(c("--tol", "0.01", "--max_iter", "30"))


# Checks
# ++++++
checks = data.table(
  case = c("dominated isoform -> 0 (default, EM run to convergence)", "dominated isoform kept > 0 with original settings (--em_tol 0.01 --em_max_iter 30)",
           "isoform-specific fragments 6:2 -> 0.75 / 0.25", "single isoform with null probability -> 1", "abundances sum to 1 per cell x gene"),
  pass = c(
    converged[["GD-B"]] < 0.005,
    original[["GD-B"]] > 0.02,
    abs(converged[["GS-A"]] - 0.75) < 0.002 & abs(converged[["GS-B"]] - 0.25) < 0.002,
    converged[["GZ-A"]] == 1,
    all(abs(tapply(converged, sub("-.*", "", names(converged)), sum) - 1) < 0.002)
  ))
checks[, result := ifelse(pass, "PASS", "FAIL")]
print(checks[, .(case, result)], row.names = FALSE)
cat("GD-B: converged", converged[["GD-B"]], "| original settings", original[["GD-B"]], "\n")
cat(sum(checks$result == "PASS"), "/", nrow(checks), "checks passed\n")
fwrite(checks[, .(case, result)], file.path(OUT, "checks.tsv"), sep = "\t")
