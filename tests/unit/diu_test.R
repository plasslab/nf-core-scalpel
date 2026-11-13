#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(Matrix, quietly = TRUE))

args = commandArgs(TRUE)
FUNCTIONS = args[1]
OUT = args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(FUNCTIONS)
set.seed(42)


# Simulated cells: isoform usage varies from cell to cell (Dirichlet), UMIs per cell ~ Poisson
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rdirichlet1 = function(a){ x = rgamma(length(a), a); x / sum(x) }
sim_gene = function(gene, p1, p2, n1, n2, conc = 2, depth = 15, depth2 = depth){
    n = n1 + n2
    pr = rbind(t(replicate(n1, rdirichlet1(conc * p1))), t(replicate(n2, rdirichlet1(conc * p2))))
    umis = c(rpois(n1, depth), rpois(n2, depth2))
    m = t(sapply(seq_len(n), function(i) rmultinom(1, umis[i], pr[i, ])))
    rownames_ = paste0(gene, "***", gene, "-", seq_along(p1), "01")
    return(Matrix(t(m), sparse = TRUE, dimnames = list(rownames_, NULL)))
}
N1 = 300
N2 = 300
null_genes = lapply(sprintf("N%03d", 1:150), function(g) sim_gene(g, c(0.5, 0.3, 0.2), c(0.5, 0.3, 0.2), N1, N2))
#- expression change only (4x more UMIs in group A), same isoform usage...
de_genes = lapply(sprintf("E%03d", 1:60), function(g) sim_gene(g, c(0.5, 0.3, 0.2), c(0.5, 0.3, 0.2), N1, N2, depth = 30, depth2 = 7))
diff_genes = lapply(sprintf("D%03d", 1:30), function(g) sim_gene(g, c(0.7, 0.3), c(0.4, 0.6), N1, N2))
#- tiny change (0.02) with a large number of UMIs: X2 significant, below --diu_min_delta...
tiny_gene = sim_gene("TINY", c(0.51, 0.49), c(0.49, 0.51), N1, N2, conc = 1000, depth = 500)
counts = do.call(rbind, c(null_genes, de_genes, diff_genes, list(tiny_gene)))
colnames(counts) = paste0("cell", seq_len(ncol(counts)))
groups = rep(c("A", "B"), c(N1, N2))

#- 3' ends: isoform -101 proximal, -201 distal (+ strand)...
ends = data.table(gene_name = sub("\\*\\*\\*.*", "", rownames(counts)), transcript_name = sub(".*\\*\\*\\*", "", rownames(counts)), strand = "+")
ends[, three_end := 1000 + 500 * (as.integer(substr(transcript_name, nchar(transcript_name) - 2, nchar(transcript_name) - 2)) - 1)]
#- terminal exon start: shared by the isoforms of the D genes (tandem 3' UTRs), distinct for the other genes...
ends[, te5 := ifelse(grepl("^D", gene_name), 900, three_end - 100)]

chisq = FindIsoforms(counts, groups, ends = ends, test = "chisq", cores = 2)
perm = FindIsoforms(counts, groups, ends = ends, test = "permutation", n_perm = 2000, cores = 3)
perm1 = FindIsoforms(counts, groups, ends = ends, test = "permutation", n_perm = 2000, cores = 1)
fpr = function(res, prefix = "^N") res[grepl(prefix, gene_name), .(p = p.value[1]), by = gene_name][, mean(p < 0.05, na.rm = TRUE)]
power = function(res) res[grepl("^D", gene_name), .(s = significant[1]), by = gene_name][, mean(s)]


# Confounding by sample: isoform usage differs between samples (batch), not between clusters
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- sample S1 (usage 0.7/0.3): 75% cluster A ; sample S2 (usage 0.4/0.6): 25% cluster A...
batch = lapply(sprintf("B%03d", 1:40), function(g) cbind(sim_gene(g, c(0.7, 0.3), c(0.7, 0.3), 200, 200), sim_gene(g, c(0.4, 0.6), c(0.4, 0.6), 200, 200)))
counts_b = do.call(rbind, batch)
colnames(counts_b) = paste0("cell", seq_len(ncol(counts_b)))
samples_b = rep(c("S1", "S2"), each = 400)
clusters_b = c(rep(c("A", "B"), c(300, 100)), rep(c("A", "B"), c(100, 300)))
perm_strat = FindIsoforms(counts_b, clusters_b, strata = samples_b, test = "permutation", n_perm = 2000, cores = 2)
perm_nostrat = FindIsoforms(counts_b, clusters_b, strata = NULL, test = "permutation", n_perm = 2000, cores = 2)
fpr_b = function(res) res[, .(p = p.value[1]), by = gene_name][, mean(p < 0.05, na.rm = TRUE)]

#- only single-isoform genes (e.g. few selected cells): no test, explicit message...
single = FindIsoforms(counts[c(1, 4), ], groups, test = "chisq")

#- replication: clusters A/B in 2 samples; gene RC changes the same way in both, gene RO in opposite directions...
cmp2 = function(g, a1, b1, a2, b2) cbind(sim_gene(g, a1, b1, 200, 200, conc = 50, depth = 20), sim_gene(g, a2, b2, 200, 200, conc = 50, depth = 20))
rep_counts = rbind(cmp2("RC", c(0.7, 0.3), c(0.3, 0.7), c(0.7, 0.3), c(0.3, 0.7)), cmp2("RO", c(0.9, 0.1), c(0.1, 0.9), c(0.35, 0.65), c(0.65, 0.35)))
colnames(rep_counts) = paste0("cell", seq_len(ncol(rep_counts)))
rep_res = FindIsoforms(rep_counts, rep(rep(c("A", "B"), each = 200), 2), strata = rep(c("S1", "S2"), each = 400), test = "chisq")
genes = function(res) unique(res[!is.na(p.value)], by = c("comparison", "gene_name"))
rep_genes = genes(rep_res)
genes_chisq = genes(chisq)

#- > 2 groups: each group vs the rest...
multi = FindIsoforms(counts[, 1:600], rep(c("A", "B", "C"), 200), test = "chisq", cores = 2)


# Checks
# ++++++
bh_gene = unique(chisq[!is.na(p.value), .(gene_name, p.value, p.value.adjusted)])
d1 = chisq[gene_name == "D001"]
checks = data.table(
  check = c(
    sprintf("chisq: false positive rate on overdispersed null genes is inflated (%.2f)", fpr(chisq)),
    sprintf("permutation: false positive rate on null genes is calibrated (%.2f <= 0.10)", fpr(perm)),
    sprintf("permutation: true usage changes detected (power %.2f >= 0.9)", power(perm)),
    sprintf("permutation: expression change without usage change not called (%.2f <= 0.10)", fpr(perm, "^E")),
    sprintf("permutation: same p-values with 1 or 3 cores"),
    "BH adjustment once per gene (one adjusted p-value per gene, = p.adjust of the gene p-values)",
    "delta_PU = PU.group1 - PU.group2 and usage sums to 1 per gene and group",
    sprintf("delta_3end < 0 when group 1 uses the proximal 3' end more (D001: %.0f nt, expected -150)", d1$delta_3end[1]),
    "tiny usage change (0.02): X2 p-value < 0.05 but not significant (--diu_min_delta 0.1)",
    sprintf("permutation stratified by sample: batch between samples not called (%.2f <= 0.10)", fpr_b(perm_strat)),
    sprintf("permutation not stratified: batch confounded with clusters called (%.2f)", fpr_b(perm_nostrat)),
    "> 2 groups: one comparison per group vs the rest (A_vs_rest, B_vs_rest, C_vs_rest)",
    "per-group count columns kept (A, B) as in the original DIU table",
    "no gene with at least 2 isoforms: explicit message instead of an error",
    "event type: tandem 3' UTR (shared terminal exon) vs alternative terminal exon",
    "direction for tandem events only (D genes: group A uses the proximal 3' end more -> shortening)",
    "dominant isoform of each group and isoform switch (D001: -101 in A, -201 in B)",
    "replication: change in the same direction in both samples -> replicated & significant",
    "replication: opposite changes in the 2 samples -> not replicated, not significant (despite p < 0.05)",
    "one row per isoform, gene-level columns identical on the rows of a gene",
    "isoforms of a gene from proximal to distal 3' end, tested flag, counts summing to the gene UMIs, sorted by adjusted p-value, gene & isoform first"),
  pass = c(
    fpr(chisq) > 0.3,
    fpr(perm) <= 0.10,
    power(perm) >= 0.9,
    fpr(perm, "^E") <= 0.10,
    isTRUE(all.equal(perm1[order(gene.tr), p.value], perm[order(gene.tr), p.value])),
    bh_gene[, uniqueN(p.value.adjusted), by = gene_name][, all(V1 == 1)] &&
      isTRUE(all.equal(bh_gene$p.value.adjusted, p.adjust(bh_gene$p.value, method = "BH"))),
    isTRUE(all.equal(chisq$delta_PU, chisq$PU.group1 - chisq$PU.group2)) &&
      chisq[, .(s1 = sum(PU.group1), s2 = sum(PU.group2)), by = gene_name][, all(abs(s1 - 1) < 1e-9 & abs(s2 - 1) < 1e-9)],
    all(abs(d1$delta_3end + 150) < 40),
    chisq[gene_name == "TINY", p.value[1] < 0.05 & !significant[1] & all(thr.var)],
    fpr_b(perm_strat) <= 0.10,
    fpr_b(perm_nostrat) > 0.5,
    setequal(unique(multi$comparison), c("A_vs_rest", "B_vs_rest", "C_vs_rest")),
    all(c("A", "B") %in% colnames(chisq)),
    identical(colnames(single), "output"),
    all(genes_chisq[grepl("^D", gene_name)]$event_type == "tandem_3UTR") & all(genes_chisq[grepl("^N", gene_name)]$event_type == "alternative_terminal_exon", na.rm = TRUE),
    all(genes_chisq[grepl("^D", gene_name) & significant == TRUE]$direction == "shortening") & all(is.na(genes_chisq[grepl("^N", gene_name)]$delta_3end)),
    genes_chisq[gene_name == "D001", top_isoform.group1 == "D001-101" & top_isoform.group2 == "D001-201" & switch],
    rep_genes[gene_name == "RC", isTRUE(replicated) & significant & n_samples == 2],
    rep_genes[gene_name == "RO", p.value < 0.05 & isFALSE(replicated) & !significant],
    chisq[, .N, by = .(comparison, gene.tr)][, all(N == 1)] & chisq[, .(n = uniqueN(paste(p.value.adjusted, significant, event_type, replicated))), by = .(comparison, gene_name)][, all(n == 1)],
    identical(chisq[gene_name == "D001"]$transcript_name, c("D001-101", "D001-201")) & all(chisq[gene_name == "D001"]$tested) &
      chisq[gene_name == "D001", isTRUE(all.equal(sum(counts.group1), umis.group1[1]))] &
      chisq[!is.na(p.value.adjusted), !is.unsorted(p.value.adjusted)] & identical(colnames(chisq)[1:2], c("gene_name", "transcript_name"))))
fwrite(checks, file.path(OUT, "checks.tsv"), sep = "\t")
print(checks)
quit(status = ifelse(all(checks$pass), 0, 1))
