# [Script for Function needed in Scalpel analysis]
# ------------------------------------------------
# Differential isoform usage (DIU) functions, from plasslab/SCALPEL src/scalpel_library.R (FindIsoforms)
#   - test "chisq" (default, original): X2 test per gene on the counts aggregated by group (isoforms x groups)
#   - test "permutation": same X2 statistic, null distribution by permutation of the cells between the
#     compared groups (within samples when clusters are compared): accounts for the cell-to-cell variability
#     of the isoform usage (cells are the units, not the UMIs).
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(Matrix, quietly = TRUE))


# Isoform 3' ends (strand-aware) from the isoform annotation (<outdir>/isoform/*.annots.gz)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Read_isoformEnds = function(annot_files){
    if(length(annot_files) == 0) return(NULL)
    annots = rbindlist(lapply(annot_files, function(f) fread(f, select = c("start", "end", "strand", "gene_name", "transcript_name", "exon_number"))))
    #- terminal exon (exon_number 1): 3' end = end (+ strand) / start (- strand)...
    annots = annots[exon_number == 1]
    annots[, three_end := ifelse(strand == "+", end, start)]
    annots = unique(annots[, .(gene_name, transcript_name, strand, three_end)], by = c("gene_name", "transcript_name"))
    #- feature names as in the Seurat object ("_" replaced by "-")...
    annots[, `:=`(gene_name = gsub("_", "-", gene_name), transcript_name = gsub("_", "-", transcript_name))]
    return(annots)
}


# Pearson X2 statistic of 2 x I tables, one table per row of S1 (counts of group 1), tot = isoform totals
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
X2_stat = function(S1, tot){
    S1 = matrix(S1, ncol = length(tot))
    N = sum(tot)
    S2 = sweep(-S1, 2, tot, "+")
    E1 = outer(rowSums(S1), tot) / N
    E2 = outer(N - rowSums(S1), tot) / N
    return(rowSums((S1 - E1)^2 / E1 + (S2 - E2)^2 / E2))
}


# Permutation p-value of the X2 statistic: cells permuted between the 2 groups (within strata)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- X: cells x isoforms counts (cells expressing the gene), lab: 1/2, strata: sample of each cell (or NULL)
#- sequential stopping (Besag & Clifford 1991): stop after `h` permuted statistics >= observed
Perm_pvalue = function(X, lab, strata = NULL, n_perm = 10000, batch = 500, h = 20){
    tot = colSums(X)
    obs = X2_stat(colSums(X[lab == 1, , drop = FALSE]), tot)
    blocks = if(is.null(strata)) list(seq_along(lab)) else split(seq_along(lab), strata)
    exceed = 0
    done = 0
    while(done < n_perm & exceed < h){
        b = min(batch, n_perm - done)
        #- membership of group 1 for each permutation (b x cells)...
        P = matrix(0, nrow = b, ncol = length(lab))
        for(i in seq_len(b)){
            perm = lab
            for(idx in blocks) perm[idx] = lab[idx][sample.int(length(idx))]
            P[i, ] = perm == 1
        }
        exceed = exceed + sum(X2_stat(P %*% X, tot) >= obs - 1e-9)
        done = done + b
    }
    return(list(statistic = obs, p.value = (exceed + 1) / (done + 1), n_perm = done))
}


# Test of one gene between 2 groups (agg: isoforms x 2 aggregated counts)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Test_gene = function(agg, X, lab, strata, test, threshold_abund, threshold_var, n_perm){
    pu = sweep(agg, 2, pmax(colSums(agg), 1e-12), "/")
    res = data.table(
        thr.abund = apply(pu, 1, function(x) all(x < threshold_abund)),
        thr.var = abs(pu[, 1] - pu[, 2]) < threshold_var,
        statistic = NA_real_, p.value = NA_real_, n_perm = NA_integer_)
    keep = which(!res$thr.abund)
    #- X2 test if more than 1 isoform following the abundance filter & X2 conditions (expected counts >= 5)...
    if(length(keep) > 1){
        numTab = agg[keep, , drop = FALSE]
        chi2_res = suppressWarnings(chisq.test(numTab))
        if(all(round(chi2_res$expected, 0) >= 5)){
            if(test == "chisq"){
                res[, `:=`(statistic = unname(chi2_res$statistic), p.value = chi2_res$p.value)]
            } else {
                Xg = X[, keep, drop = FALSE]
                expressed = rowSums(Xg) > 0
                pr = Perm_pvalue(Xg[expressed, , drop = FALSE], lab[expressed], if(is.null(strata)) NULL else strata[expressed], n_perm)
                res[, `:=`(statistic = pr$statistic, p.value = pr$p.value, n_perm = as.integer(pr$n_perm))]
            }
        }
    }
    return(res)
}


# Differential isoform usage between 2 groups of cells (group1 vs group2)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- counts: isoforms x cells (features gene***isoform), lab: 1/2 per cell, strata: sample per cell (or NULL)
Compare_groups = function(counts, lab, strata, features, ends, test, threshold_abund, threshold_var, n_perm, cores, seed){
    agg = cbind(Matrix::rowSums(counts[, lab == 1, drop = FALSE]), Matrix::rowSums(counts[, lab == 2, drop = FALSE]))
    genes = split(seq_len(nrow(features)), features$gene_name)
    genes = genes[lengths(genes) > 1]
    if(length(genes) == 0) return(NULL)
    tcounts = if(test == "permutation") Matrix::t(counts) else NULL

    #- one test per gene (seed per gene: reproducible whatever the number of cores)...
    res = rbindlist(parallel::mclapply(seq_along(genes), function(i){
        set.seed(seed + i)
        idx = genes[[i]]
        X = if(test == "permutation") as.matrix(tcounts[, idx, drop = FALSE]) else NULL
        cbind(features[idx], Test_gene(agg[idx, , drop = FALSE], X, lab, strata, test, threshold_abund, threshold_var, n_perm))
    }, mc.cores = cores))

    #- isoform usage (PU) per group, change of usage & weighted 3' end shift (+: group1 uses more distal 3' ends)...
    res = cbind(res, counts.group1 = agg[res$row, 1], counts.group2 = agg[res$row, 2])
    res[, `:=`(PU.group1 = counts.group1 / pmax(sum(counts.group1), 1e-12), PU.group2 = counts.group2 / pmax(sum(counts.group2), 1e-12)), by = gene_name]
    res[, delta_PU := PU.group1 - PU.group2]
    if(!is.null(ends)){
        res = merge(res, ends, by = c("gene_name", "transcript_name"), all.x = TRUE, sort = FALSE)
        res[, dist_3end := if(all(is.na(three_end))) NA_real_ else as.numeric(ifelse(strand == "+", three_end - min(three_end, na.rm = TRUE), max(three_end, na.rm = TRUE) - three_end)), by = gene_name]
        res[, `:=`(end_rank = frank(dist_3end, ties.method = "dense", na.last = "keep"),
                   delta_3end = sum(delta_PU * dist_3end)), by = gene_name]
    }

    #- BH adjustment once per gene (one test per gene)...
    padj = unique(res[!is.na(p.value), .(gene_name, p.value)])
    padj[, p.value.adjusted := p.adjust(p.value, method = "BH")]
    res = merge(res, padj[, .(gene_name, p.value.adjusted)], by = "gene_name", all.x = TRUE, sort = FALSE)
    res[, max_delta_PU := max(abs(delta_PU[!thr.abund]), 0), by = gene_name]
    return(res)
}


# FindIsoforms: all comparisons (2 groups: group1 vs group2; > 2 groups: each group vs the rest)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FindIsoforms = function(counts, groups, strata = NULL, ends = NULL, test = "chisq",
                        threshold_pval = 0.05, threshold_abund = 0.05, threshold_var = 0.1,
                        n_perm = 10000, cores = 1, seed = 1){
    if(!(test %in% c("chisq", "permutation"))) stop("Unknown DIU test '", test, "' (chisq or permutation)")
    groups = as.character(groups)
    lv = sort(unique(groups))
    if(length(lv) < 2) stop("Provide more conditions (levels > 1) for differential isoform usage analysis")
    features = data.table(row = seq_len(nrow(counts)), gene.tr = rownames(counts))
    features[, c("gene_name", "transcript_name") := tstrsplit(gene.tr, "***", fixed = TRUE)]
    comparisons = if(length(lv) == 2) list(c(lv[1], lv[2])) else lapply(lv, function(k) c(k, "rest"))

    all_results = lapply(comparisons, function(cmp){
        message("Comparison ", cmp[1], " vs ", cmp[2], " (", test, ")...")
        lab = ifelse(groups == cmp[1], 1, ifelse(cmp[2] == "rest" | groups == cmp[2], 2, NA))
        cells = which(!is.na(lab))
        #- strata (samples) used only if a sample holds cells of both groups (else not stratified)...
        st = if(is.null(strata)) NULL else as.character(strata)[cells]
        if(!is.null(st)){
            both = tapply(lab[cells], st, function(x) length(unique(x)) == 2)
            if(!any(both)){
                message("  groups not shared by any sample: permutation not stratified by sample")
                st = NULL
            }
        }
        res = Compare_groups(counts[, cells, drop = FALSE], lab[cells], st, features, ends, test,
                             threshold_abund, threshold_var, n_perm, cores, seed)
        if(is.null(res)) return(NULL)
        res[, `:=`(comparison = paste0(cmp[1], "_vs_", cmp[2]), group1 = cmp[1], group2 = cmp[2])]
        return(res)
    })
    all_results = rbindlist(all_results, fill = TRUE)
    if(nrow(all_results) == 0){
        message("No gene with at least 2 isoforms: no differential isoform usage test")
        return(data.table(output = "No gene with at least 2 isoforms for differential isoform usage analysis"))
    }

    #- counts aggregated by group (one column per group, as the original DIU table)...
    agg_all = sapply(lv, function(k) Matrix::rowSums(counts[, groups == k, drop = FALSE]))
    agg_all = data.table(gene.tr = rownames(counts), matrix(agg_all, ncol = length(lv), dimnames = list(NULL, lv)))
    all_results = merge(agg_all, all_results, by = "gene.tr", sort = FALSE)

    #- significant genes: adjusted p-value & change of usage of at least one tested isoform...
    all_results[, significant := !is.na(p.value.adjusted) & p.value.adjusted < threshold_pval & max_delta_PU >= threshold_var]
    all_results[, test := test]
    setorder(all_results, comparison, p.value.adjusted, gene_name, transcript_name, na.last = TRUE)
    all_results[, row := NULL]
    first = c("comparison", "group1", "group2", "gene.tr", "gene_name", "transcript_name", lv,
              "counts.group1", "counts.group2", "PU.group1", "PU.group2", "delta_PU", "strand", "three_end", "dist_3end", "end_rank",
              "thr.abund", "thr.var", "statistic", "p.value", "p.value.adjusted", "max_delta_PU", "delta_3end", "significant", "test", "n_perm")
    setcolorder(all_results, intersect(first, colnames(all_results)))
    return(all_results)
}
