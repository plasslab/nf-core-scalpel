
suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(stringr, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
options(width = 300)


parser <- ArgumentParser(description = "Calculate fragment probability for each read and isoform combination")
parser$add_argument("--reads", required = TRUE, help = "Path to reads file")
parser$add_argument("--probs", required = TRUE, help = "Path to probability distribution file")
parser$add_argument("--output", required = TRUE, help = "Path to output file")
parser$add_argument("--threads", type = "integer", help = "Number of threads")
args <- parser$parse_args()


#0. Opening
reads = fst::read_fst(args$reads, as.data.table = TRUE)
probs = fread(
    args$probs, 
    col.names = c("start.rdR","counts","probs_bin"),
    nThread = args$threads
)

reads %>% print()
probs %>% print()

#1. Joining
reads = left_join(reads, probs, by="start.rdR") %>%
    mutate(
        probs_bin = ifelse(is.na(probs_bin), 0, probs_bin),
        counts = ifelse(is.na(counts), 0, counts)
    )

#2. Distinct and separate tags
reads = reads %>% 
    dplyr::distinct(
        read_start,
        read_end,
        strand,
        tags,
        gene_name,
        transcript_name,
        TPM_perc,
        probs_bin
    ) %>%
    tidyr::separate(
        col="tags", 
        into=c("bc","umi"), 
        sep="::", 
        remove=T) %>%
    data.table()

#3. Grouping
reads = reads %>% 
    group_by(bc,umi,gene_name,transcript_name) %>% 
    dplyr::mutate(frag_probs_weighted = TPM_perc * (prod(probs_bin))) %>% 
    dplyr::distinct(bc,gene_name,transcript_name,umi,frag_probs_weighted,TPM_perc) %>% 
    arrange(bc,gene_name,transcript_name,umi) %>% 
    dplyr::distinct(bc,gene_name,transcript_name,umi,frag_probs_weighted) %>% 
    data.table()

#delete sequencing tags
reads$bc = str_replace(reads$bc, pattern = "XC:Z:","")
reads$bc = str_replace(reads$bc, pattern = "CB:Z:","")
reads$umi = str_replace(reads$umi, pattern = "XM:Z:","")
reads$umi = str_replace(reads$umi, pattern = "UB:Z:","")

# write output
fst::write_fst(reads, args$output)
