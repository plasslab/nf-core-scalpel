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
    #- terminal exon (exon_number 1): 3' end = end (+ strand) / start (- strand), 5' boundary te5 = start (+) / end (-)...
    annots = annots[exon_number == 1]
    annots[, `:=`(three_end = ifelse(strand == "+", end, start), te5 = ifelse(strand == "+", start, end))]
    annots = unique(annots[, .(gene_name, transcript_name, strand, three_end, te5)], by = c("gene_name", "transcript_name"))
    #- feature names as in the Seurat object ("_" replaced by "-")...
    annots[, `:=`(gene_name = gsub("_", "-", gene_name), transcript_name = gsub("_", "-", transcript_name))]
    return(annots)
}


# Permutation p-values of the X2 statistic for all tested genes at once
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- cell labels permuted between the 2 groups (within strata), the same permutations for all genes:
#  permuted counts of group 1 of all isoforms = one sparse product per batch of permutations
#- sequential stopping per gene (Besag & Clifford 1991): a gene stops after `h` permuted statistics >= observed
#- counts: tested isoforms x cells, gene: gene of each tested isoform, lab: 1/2 per cell, strata: sample per cell (or NULL)
X2_genes = function(S1, tot, gid, M){
    #- S1: permutations x isoforms counts of group 1 -> permutations x genes X2 statistics...
    Nf = as.vector(tot %*% M)[gid]
    r1 = as.matrix(S1 %*% M)[, gid, drop = FALSE]
    E1 = sweep(r1, 2, tot / Nf, "*")
    E2 = sweep(sweep(-r1, 2, Nf, "+"), 2, tot / Nf, "*")
    S2 = sweep(-S1, 2, tot, "+")
    terms = (S1 - E1)^2 / E1 + (S2 - E2)^2 / E2
    terms[!is.finite(terms)] = 0
    return(as.matrix(terms %*% M))
}

Perm_genes = function(counts, gene, lab, strata = NULL, n_perm = 10000, batch = 500, h = 20, seed = 1){
    gid = as.integer(factor(gene, levels = unique(gene)))
    G = max(gid)
    tcounts = Matrix::t(counts)
    tot = Matrix::colSums(tcounts)
    M = Matrix::sparseMatrix(i = seq_along(gid), j = gid, x = 1, dims = c(length(gid), G))
    obs = X2_genes(matrix(Matrix::colSums(tcounts[lab == 1, , drop = FALSE]), nrow = 1), tot, gid, M)[1, ]
    blocks = if(is.null(strata)) list(seq_along(lab)) else split(seq_along(lab), strata)
    exceed = numeric(G)
    done = numeric(G)
    active = seq_len(G)
    set.seed(seed)
    while(length(active) > 0){
        b = min(batch, n_perm - min(done[active]))
        #- group 1 membership of the cells for each permutation (cells x b)...
        L = matrix(0, nrow = length(lab), ncol = b)
        for(i in seq_len(b)){
            perm = lab
            for(idx in blocks) perm[idx] = lab[idx][sample.int(length(idx))]
            L[, i] = perm == 1
        }
        #- permuted X2 of the active genes...
        fa = which(gid %in% active)
        ga = match(gid[fa], active)
        S1 = t(as.matrix(Matrix::crossprod(tcounts[, fa, drop = FALSE], L)))
        x2 = X2_genes(S1, tot[fa], ga, M[fa, active, drop = FALSE])
        exceed[active] = exceed[active] + colSums(sweep(x2, 2, obs[active] * (1 - 1e-9), ">="))
        done[active] = done[active] + b
        active = active[done[active] < n_perm & exceed[active] < h]
    }
    return(data.table(gene_name = unique(gene), statistic = obs, p.value = (exceed + 1) / (done + 1), n_perm = as.integer(done)))
}


# Filters & X2 test of one gene between 2 groups (agg: isoforms x 2 aggregated counts)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Test_gene = function(agg, threshold_abund, threshold_var){
    pu = sweep(agg, 2, pmax(colSums(agg), 1e-12), "/")
    res = data.table(
        thr.abund = apply(pu, 1, function(x) all(x < threshold_abund)),
        thr.var = abs(pu[, 1] - pu[, 2]) < threshold_var,
        tested = FALSE, statistic = NA_real_, p.value = NA_real_)
    keep = which(!res$thr.abund)
    #- X2 test if more than 1 isoform following the abundance filter & X2 conditions (expected counts >= 5)...
    if(length(keep) > 1){
        chi2_res = suppressWarnings(chisq.test(agg[keep, , drop = FALSE]))
        if(all(round(chi2_res$expected, 0) >= 5)){
            res[keep, tested := TRUE]
            res[, `:=`(statistic = unname(chi2_res$statistic), p.value = chi2_res$p.value)]
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

    #- filters & X2 test per gene...
    res = rbindlist(lapply(genes, function(idx) cbind(features[idx], Test_gene(agg[idx, , drop = FALSE], threshold_abund, threshold_var))))
    res[, n_perm := NA_integer_]

    #- permutation p-values of the tested genes (genes split on the cores, same permutations for all genes)...
    if(test == "permutation" & any(res$tested)){
        tested = res[tested == TRUE, .(row, gene_name)]
        chunks = split(unique(tested$gene_name), rep_len(seq_len(cores), uniqueN(tested$gene_name)))
        perm = rbindlist(parallel::mclapply(chunks, function(g){
            sub = tested[gene_name %in% g]
            Perm_genes(counts[sub$row, , drop = FALSE], sub$gene_name, lab, strata, n_perm, seed = seed)
        }, mc.cores = cores))
        res[, c("statistic", "p.value") := NULL]
        res = merge(res, perm, by = "gene_name", all.x = TRUE, sort = FALSE, suffixes = c(".x", ""))
        res[, n_perm.x := NULL]
    }
    res[, tested := NULL]

    #- isoform usage (PU) per group, change of usage & weighted 3' end shift (+: group1 uses more distal 3' ends)...
    res = cbind(res, counts.group1 = agg[res$row, 1], counts.group2 = agg[res$row, 2])
    res[, `:=`(PU.group1 = counts.group1 / pmax(sum(counts.group1), 1e-12), PU.group2 = counts.group2 / pmax(sum(counts.group2), 1e-12)), by = gene_name]
    res[, delta_PU := PU.group1 - PU.group2]
    res[, `:=`(umis.group1 = sum(counts.group1), umis.group2 = sum(counts.group2), tested = !is.na(p.value) & !thr.abund), by = gene_name]
    #- dominant isoform of each group & isoform switch...
    res[, `:=`(top_isoform.group1 = transcript_name[which.max(PU.group1)], top_isoform.group2 = transcript_name[which.max(PU.group2)]), by = gene_name]
    res[, switch := top_isoform.group1 != top_isoform.group2]

    #- event type: tandem 3' UTR (tested isoforms share their terminal exon start: alternative polyadenylation in the
    #  same last exon) / alternative terminal exon; 3' UTR shift & direction for tandem events only...
    if(is.null(ends)) ends = data.table(gene_name = character(), transcript_name = character(), strand = character(), three_end = numeric(), te5 = numeric())
    res = merge(res, ends, by = c("gene_name", "transcript_name"), all.x = TRUE, sort = FALSE)
    res[, dist_3end := if(all(is.na(three_end))) NA_real_ else as.numeric(ifelse(strand == "+", three_end - min(three_end, na.rm = TRUE), max(three_end, na.rm = TRUE) - three_end)), by = gene_name]
    res[, end_rank := frank(dist_3end, ties.method = "dense", na.last = "keep"), by = gene_name]
    res[, event_type := {
        t5 = te5[!thr.abund]
        if(length(t5) < 2 | any(is.na(t5))) NA_character_ else if(uniqueN(t5) == 1) "tandem_3UTR" else "alternative_terminal_exon"
    }, by = gene_name]
    res[, delta_3end := if(isTRUE(event_type[1] == "tandem_3UTR")) sum((delta_PU * dist_3end)[!thr.abund]) else NA_real_, by = gene_name]
    res[, direction := fifelse(event_type == "tandem_3UTR", fifelse(delta_3end < 0, "shortening", "lengthening"), NA_character_)]

    #- BH adjustment once per gene (one test per gene)...
    padj = unique(res[!is.na(p.value), .(gene_name, p.value)])
    padj[, p.value.adjusted := p.adjust(p.value, method = "BH")]
    res = merge(res, padj[, .(gene_name, p.value.adjusted)], by = "gene_name", all.x = TRUE, sort = FALSE)
    res[, max_delta_PU := max(abs(delta_PU[!thr.abund]), 0), by = gene_name]
    return(res)
}


# Replication across samples: change of the most changing isoform of each gene in each sample
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#- a sample is evaluable for a gene if the gene has >= min_umis UMIs in both groups of the sample, and consistent if
#  the change has the same sign as the overall change and is >= threshold_var / 2
Replicate_effects = function(counts, lab, samples, res, threshold_var, min_umis = 20){
    tested = res[!is.na(p.value)]
    if(nrow(tested) == 0) return(NULL)
    #- most changing isoform (ties: the increased one)...
    top = tested[thr.abund == FALSE][order(gene_name, -abs(round(delta_PU, 9)), -delta_PU)][, .SD[1], by = gene_name, .SDcols = c("row", "transcript_name", "delta_PU")]
    gene_rows = tested[, .(row), by = gene_name]
    eff = rbindlist(lapply(sort(unique(samples)), function(s){
        c1 = which(lab == 1 & samples == s)
        c2 = which(lab == 2 & samples == s)
        if(length(c1) == 0 | length(c2) == 0) return(NULL)
        a1 = Matrix::rowSums(counts[, c1, drop = FALSE])
        a2 = Matrix::rowSums(counts[, c2, drop = FALSE])
        tot = gene_rows[, .(t1 = sum(a1[row]), t2 = sum(a2[row])), by = gene_name]
        x = merge(top, tot, by = "gene_name")
        x[, `:=`(sample = s, d = a1[row] / pmax(t1, 1e-12) - a2[row] / pmax(t2, 1e-12), ok = t1 >= min_umis & t2 >= min_umis)]
        return(x)
    }))
    if(nrow(eff) == 0) return(NULL)
    eff[, consistent := ok & sign(d) == sign(delta_PU) & abs(d) >= threshold_var / 2]
    eff = eff[, .(replication_isoform = transcript_name[1], n_samples = sum(ok), n_samples_consistent = sum(consistent),
                  delta_PU_by_sample = paste0(sample, ":", ifelse(ok, round(d, 3), "NA"), collapse = ";")), by = gene_name]
    eff[, replicated := fifelse(n_samples >= 2, n_samples_consistent == n_samples, NA)]
    return(eff)
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
        #- replication of the change in each sample (clusters compared across samples)...
        rep = if(is.null(strata)) NULL else Replicate_effects(counts[, cells, drop = FALSE], lab[cells], as.character(strata)[cells], res, threshold_var)
        if(is.null(rep)) rep = data.table(gene_name = character(), replication_isoform = character(), n_samples = integer(), n_samples_consistent = integer(), delta_PU_by_sample = character(), replicated = logical())
        res = merge(res, rep, by = "gene_name", all.x = TRUE, sort = FALSE)
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
    #- ... and, with clusters compared across samples, a change in the same direction in each evaluable sample...
    all_results[, significant := !is.na(p.value.adjusted) & p.value.adjusted < threshold_pval & max_delta_PU >= threshold_var & (is.na(replicated) | replicated)]
    all_results[, test := test]
    all_results[, row := NULL]
    #- one row per isoform: gene-level columns first (repeated), then the isoform columns...
    first = c("gene_name", "transcript_name", "comparison", "group1", "group2", "event_type", "direction", "delta_3end",
              "tested", "counts.group1", "counts.group2", "PU.group1", "PU.group2", "delta_PU",
              "p.value.adjusted", "p.value", "statistic", "max_delta_PU", "replicated", "replication_isoform",
              "n_samples", "n_samples_consistent", "delta_PU_by_sample", "switch", "top_isoform.group1", "top_isoform.group2",
              "umis.group1", "umis.group2", "dist_3end", "end_rank", "strand", "three_end", "te5", "thr.abund", "thr.var",
              "gene.tr", lv, "test", "n_perm", "significant")
    setcolorder(all_results, intersect(first, colnames(all_results)))
    #- genes by adjusted p-value & usage change, isoforms of a gene from the most proximal to the most distal 3' end...
    all_results[, o_d := -max_delta_PU]
    setorder(all_results, comparison, p.value.adjusted, o_d, gene_name, dist_3end, transcript_name, na.last = TRUE)
    all_results[, o_d := NULL]
    return(all_results)
}
