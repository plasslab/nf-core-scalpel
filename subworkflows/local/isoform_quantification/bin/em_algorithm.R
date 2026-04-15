
suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))


parser <- ArgumentParser(description = "EM algorithm for isoform quantification")
parser = ArgumentParser(description='Probability of unique transcript')
parser$add_argument('--cell', type="character", help='path of cell file')
parser$add_argument('--output_path', type="character", help='output path')
args = parser$parse_args()


#Define EM algorithm function
#============================

em_algorithm <- function(bc_gene, tx_id, umi_id, prob, transcript_name, max_iter = 30L, tol = 1e-2) {

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

#- arrange ordering...
setorder(df, bc_gene, umi_local, tx_id)

# avoid null isoform bulk probability values, which can cause issues with the EM algorithm. 
# Adding a small constant to ensure all probabilities are non-zero.
df$frag_probs_weighted = df$frag_probs_weighted + 1e-9

#- perfom EM algorithm for each cell-gene combination...
cell = df[, em_algorithm(bc_gene, tx_id, umi_local, frag_probs_weighted, transcript_name), by=.(bc_gene, bc, gene_name)]

#- writing output...
write_fst(cell, args$output_path)
