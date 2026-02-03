#!/usr/bin/env Rscript

require(argparse)
require(dplyr)

# Parser setup
parser <- argparse::ArgumentParser(description = "Collapse transcripts by isoform end profile")
parser$add_argument("--input", "-i", required = TRUE, help = "Input BED file (preprocessed GTF)")
parser$add_argument("--output", "-o", required = TRUE, help = "Output file with collapsed_transcripts column")
parser$add_argument("--distance_3end", "-d3", type = "integer", default = 30, 
                    help = "Maximum distance between 3' ends (default: 30bp)")
parser$add_argument("--distance_profile", "-dp", type = "integer", default = 600,
                    help = "Maximum relative transcriptomic distance (default: 600bp)")
args <- parser$parse_args()

# Load the function
groupIsoformEndProfile <- function(gtf_tab_with_rel_coords, distance_3end = 30, distance_profile = 600) {
    cluster_by_proximity <- function(x, distance) {
        if (length(x) == 0) return(integer(0))
        if (length(x) == 1) return(1L)
        idx <- order(x)
        x_sorted <- x[idx]
        res <- base::cumsum(c(1L, diff(x_sorted) > distance))
        res_out <- integer(length(x))
        res_out[idx] <- res
        res_out
    }
    
    last_exons <- gtf_tab_with_rel_coords %>%
        dplyr::filter(startR == 0) %>%
        dplyr::distinct(gene_id, transcript_name, .keep_all = TRUE) %>%
        dplyr::mutate(three_prime_end = ifelse(strand == "+", end, start)) %>%
        dplyr::group_by(gene_id, strand) %>%
        dplyr::mutate(cluster_3end = cluster_by_proximity(three_prime_end, distance_3end)) %>%
        dplyr::group_by(gene_id, strand, cluster_3end) %>%
        dplyr::mutate(
            cluster_3end_ref = ifelse(strand == "+", max(three_prime_end), min(three_prime_end)),
            offset_to_ref = abs(cluster_3end_ref - three_prime_end)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(gene_id, transcript_name, cluster_3end, three_prime_end, cluster_3end_ref, offset_to_ref)
    
    exons_with_rel_boundary <- gtf_tab_with_rel_coords %>%
        dplyr::left_join(last_exons, by = c("gene_id", "transcript_name"), relationship = "many-to-one") %>%
        dplyr::filter(!is.na(cluster_3end)) %>%
        dplyr::mutate(
            startR_adj = startR + offset_to_ref,
            endR_adj = endR + offset_to_ref
        ) %>%
        dplyr::filter(startR_adj < distance_profile) %>%
        dplyr::mutate(
            rel_5prime_boundary = endR_adj + 1,
            rel_5prime_boundary_capped = pmin(rel_5prime_boundary, distance_profile)
        )
    
    profile_clusters <- exons_with_rel_boundary %>%
        dplyr::group_by(gene_id, cluster_3end, transcript_name) %>%
        dplyr::arrange(exon_number) %>%
        dplyr::summarise(
            profile_signature = paste(rel_5prime_boundary_capped, collapse = "_"),
            n_exons_in_window = n(),
            strand = dplyr::first(strand),
            .groups = "drop"
        ) %>%
        dplyr::group_by(gene_id, cluster_3end, profile_signature) %>%
        dplyr::mutate(
            profile_group = dplyr::cur_group_id(),
            n_isoforms_in_profile = n()
        ) %>%
        dplyr::ungroup()
    
    collapsed_lookup <- profile_clusters %>%
        dplyr::group_by(profile_group) %>%
        dplyr::mutate(transcript_collapseds = paste0(sort(unique(transcript_name)), collapse = ",")) %>%
        dplyr::ungroup() %>%
        dplyr::select(transcript_name, transcript_collapseds, n_isoforms_in_profile) %>%
        dplyr::distinct()
    
    result <- gtf_tab_with_rel_coords %>%
        dplyr::left_join(collapsed_lookup, by = "transcript_name", relationship = "many-to-one")
    
    return(result)
}

# Read input file
base::cat("Reading input file:", args$input, "\n")
gtf_data <- utils::read.csv(args$input, sep = "\t", header = TRUE)

# Apply function
base::cat("Processing transcripts...\n")
result <- groupIsoformEndProfile(
    gtf_data, 
    distance_3end = args$distance_3end,
    distance_profile = args$distance_profile
)

nb.collapseds = distinct(result, transcript_collapseds, n_isoforms_in_profile) %>%
    dplyr::filter(n_isoforms_in_profile>1) %>%
    summarise(nb.collapseds = sum(n_isoforms_in_profile))

check.stats1 = distinct(result, gene_name, transcript_collapseds, n_isoforms_in_profile) %>% 
    as_tibble() %>%
    summarise(
        nb.genes = n_distinct(gene_name), 
        nb.isoforms = sum(n_isoforms_in_profile),
        nb.collapseds = nb.collapseds$nb.collapseds
    )

stats.dat = paste0(result$seqnames[1], "\t", check.stats1$nb.genes, "\t", check.stats1$nb.isoforms, "\t", check.stats1$nb.collapseds)
system(
    paste0(
        "echo ", 
        stats.dat, 
        " > stats.txt")
)
       
    
# stats_output <- c(
#     # paste0("- Statistic about annotation processing - ", result$seqnames[1]),
#     paste("Number of genes processed: ", check.stats1$nb.genes),
#     paste("Number of isoforms before collapsing: ", check.stats1$nb.isoforms),
#     paste("Number of collapsed isoforms: ", check.stats1$nb.collapseds)
# )
# writeLines(stats_output, con = "stats.txt")

# Write output
base::cat("Writing output to:", args$output, "\n")
utils::write.table(result, file = args$output, row.names = FALSE, quote = FALSE, sep = "\t")
base::cat("Done!\n")
