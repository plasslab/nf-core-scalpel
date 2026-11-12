
suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))


parser <- ArgumentParser(description = "EM algorithm for isoform quantification")
parser = ArgumentParser(description='Probability of unique transcript')
parser$add_argument('--cell', type="character", help='path of cell file')
parser$add_argument('--output_path', type="character", help='output path')
parser$add_argument('--max_iter', type="character", default="1000", help='maximum number of EM iterations (original SCALPEL: 30)')
parser$add_argument('--tol', type="character", default="1e-6", help='stop when no abundance changes by more than tol (original SCALPEL: 0.01)')
args = parser$parse_args()

MAX_ITER = as.integer(args$max_iter)
TOL = as.numeric(args$tol)


#Define EM algorithm function
#============================

em_algorithm <- function(bc_gene, tx_id, umi_id, prob, transcript_name, max_iter = MAX_ITER, tol = TOL) {

    # remap to local contiguous indices within this group
    tx_unique = unique(tx_id)
    tx_local  = match(tx_id, tx_unique)
    umi_local = match(umi_id, unique(umi_id))

    # transcript names in same order as tx_local indices
    tr_names = transcript_name[match(tx_unique, tx_id)]

    nb_transcripts = length(tx_unique)
    nb_umis = max(umi_local)

    # case 1: only one transcript
    if (nb_transcripts == 1L) {
        return(data.table(transcript_name = tr_names, rel_abund = 1.0))
    }

    theta = rep(1 / nb_transcripts, nb_transcripts)

    for (i in seq_len(max_iter)) {
        prior_theta = theta

        # E-step
        weighted = prob * theta[tx_local]
        denom = rowsum(weighted, umi_local, reorder = FALSE)[, 1]
        posterior = weighted / denom[umi_local]
        # M-step
        theta = rowsum(posterior, tx_local, reorder = FALSE)[, 1] / nb_umis
        # convergence
        if (max(abs(theta - prior_theta)) < tol) {
            break
        }
    }
    return(data.table(transcript_name = tr_names, rel_abund = theta %>% round(3)))
}

#- reading...
df = read_fst(args$cell, as.data.table = TRUE)

df[, bc_gene := .GRP, by = .(bc, gene_name)]
df[, tx_id := .GRP, by = .(bc, gene_name, transcript_name)]
df[, umi_local := .GRP, by = .(bc, gene_name, umi)]

# avoid null isoform bulk probability values, which can cause issues with the EM algorithm. 
# Adding a small constant to ensure all probabilities are non-zero.
#- (fragments with a null probability for all isoforms, i.e. reads beyond the learned distances, are kept:
#-  compatible with a single isoform, they remain a molecule of this isoform)
df$frag_probs_weighted = df$frag_probs_weighted + 1e-9

#- arrange ordering...
setorder(df, bc_gene, umi_local, tx_id)

#- perfom EM algorithm for each cell-gene combination...
cell = df[, em_algorithm(bc_gene, tx_id, umi_local, frag_probs_weighted, transcript_name), by=.(bc_gene, bc, gene_name)]
if(nrow(cell) == 0){
    cell = data.table(bc_gene = integer(0), bc = character(0), gene_name = character(0), transcript_name = character(0), rel_abund = numeric(0))
}

#- writing output...
write_fst(cell, args$output_path)
