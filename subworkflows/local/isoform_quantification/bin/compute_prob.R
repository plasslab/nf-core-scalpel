#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse, quietly =TRUE))
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(scales, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(fst, quietly = TRUE))
options(width = 500)

parser <- ArgumentParser(description = "Merge reads with genomic annotations")
parser$add_argument("--repo", required = TRUE, help = "Path to reads repository")
parser$add_argument("--gene_frac", required = TRUE, help = "Fraction of expressed gene for filtering")
parser$add_argument("--bins", required = TRUE, help = "Bins size for probability distribution")
parser$add_argument("--threads", type = "integer", help = "Number of threads")
parser$add_argument("--OUTPUT_PROB", required = TRUE, help = "Path to output probability file")
parser$add_argument("--OUTPUT_PDF", required = TRUE, help = "Path to output PDF file")
args <- parser$parse_args()

BINS = args$bins %>% as.numeric()
THRESHOLD_FRAC = args$gene_frac %>% as.character()


# Read files
message("Reading files...")
reads = lapply(
    list.files(
        args$repo, 
        full.names = TRUE, 
        pattern = "*.ip_filtered.fst"
    ), function(input_file) {
        message("Reading file: ", input_file)
        fst::read_fst(input_file) %>%
            dplyr::group_by(gene_name) %>%
            dplyr::filter(n_distinct(transcript_name)==1) %>%
            data.table() %>%
            return(.)
    }
) %>% 
rbindlist()

gene_counts = reads %>%
    distinct(seqnames, read_start, read_end, strand, tags, gene_name, transcript_name) %>%
    group_by(gene_name) %>%
    summarise(nb.reads = n()) %>%
    data.table()
gene_quantiles = stats::quantile(gene_counts$nb.reads, seq(0,1,0.01))
gene_tokeep = dplyr::filter(gene_counts, nb.reads < gene_quantiles[[THRESHOLD_FRAC]])$gene_name

reads = reads %>%
    filter(gene_name %in% gene_tokeep) %>%
    data.table()

#- get the number of distinct reads at each 3'end position
read_tab = distinct(reads, start.rdR, read_start, read_end) %>%
  group_by(start.rdR) %>%
  reframe(read.counts = n()) %>%
  arrange(start.rdR) %>%
  stats::na.omit() %>%
  data.table()
colnames(read_tab) = c("transcriptomic_distance", "counts")

#let's set interval axis
part_neg = c(
  seq(0,read_tab$transcriptomic_distance[1],-BINS),
  read_tab$transcriptomic_distance[1]
) %>%
unique() %>% 
rev()
part_pos = c(
  seq(0, max(read_tab$transcriptomic_distance), BINS), 
  max(read_tab$transcriptomic_distance)
) %>% 
unique()
intervals = unique(c(part_neg,part_pos))

intervals_counts = list()
intervals_idx = list()
intervals_probs = list()
for(i in 1:(length(intervals)-1)){
  left_born = intervals[i]
  right_born = intervals[i+1]
  #index
  intervals_idx[[i]] = seq(left_born, right_born)
  #counts
  intervals_counts[[i]] = rep(
    (read_tab %>% 
      dplyr::filter(transcriptomic_distance>= left_born & 
        transcriptomic_distance < right_born
      )
    )$counts %>% 
    sum(), length(intervals_idx[[i]])
  )
}
interval.tab = data.table(
  transcriptomic_distance = unlist(intervals_idx), 
  counts_cum = unlist(intervals_counts)
) %>% 
distinct(transcriptomic_distance,.keep_all = T)
interval.tab$probs_bin = (interval.tab$counts_cum / sum(interval.tab$counts_cum)) * 100
interval.tab = left_join(interval.tab, read_tab)
interval.tab$counts = interval.tab$counts %>% tidyr::replace_na(0)

#let's normalize for visualization
interval.tab$probability_scaled = interval.tab$probs_bin %>% 
  scales::rescale(to = c(0,max(interval.tab$counts)))

# [Writing]
ggplot(interval.tab) +
  geom_area(aes(transcriptomic_distance,counts), fill="cornflowerblue", size=0.5) +
  geom_line(aes(transcriptomic_distance,probability_scaled), color="red",size=1) +
  theme_classic() +
  ggtitle("Reads distribution on transcriptomic space")
  
ggsave(
  args$OUTPUT_PDF, 
  scale = 1, 
  device = "pdf", 
  units = "in", 
  width = 11.69, 
  height = 8.27
)
dev.off()

fwrite(
  interval.tab[,c('transcriptomic_distance','counts','probs_bin')], 
  file = args$OUTPUT_PROB, 
  sep="\t", 
  col.names = F, 
  nThread = as.numeric(args$threads)
)


