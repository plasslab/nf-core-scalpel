
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(fst)))


# Argument Parser
# +++++++++++++++
parser = ArgumentParser(description='Filtering of internal priming reads')
parser$add_argument('--reads', type='character', help='path of bed_gtf file')
parser$add_argument('--ip_reference', type="character", help='path of internal priming ref file')
parser$add_argument('--output_path', type='character', help='path of output reads file')
parser$add_argument("--threads", help = "Number of threads")
parser$add_argument('--threshold_dist', type="character", help='distance of ip position from isoform end')
args = parser$parse_args()

print(args$threads)

# Load data
# +++++++++
reads = fst::read_fst(path = args$reads, as.data.table = T)
ip_pos = fread(
  file = args$ip_reference, 
  col.names = c("seqnames", "start_ip", "end_ip", "V4", "V5", "strand_ip"),
  nThread = as.numeric(args$threads)
) %>%
  filter(seqnames %in% unique(reads$seqnames))



# Filtering of internal priming reads
# +++++++++++++++++++++++++++++++++++
if(nrow(ip_pos) == 0){
  message("No internal priming positions found in the provided annotation file.")

} else {

    reads = reads %>%
        group_by(tags, transcript_name) %>%
        mutate(fidr = cur_group_id()) %>%
        data.table()

    gr_reads = GRanges(
        seqnames = reads$seqnames,
        ranges = IRanges(start = reads$read_start, end = reads$read_end),
        strand = reads$strand,
        anno_start = reads$anno_start,
        anno_end = reads$anno_end,
        startR = reads$startR,
        endR = reads$endR,
        gene_name = reads$gene_name,
        transcript_name = reads$transcript_name,
        tags = reads$tags,
        fidr = reads$fidr
    )
    gr_ip = GRanges(
        seqnames = ip_pos$seqnames,
        ranges = IRanges(start = ip_pos$start_ip, end = ip_pos$end_ip),
        strand = ip_pos$strand_ip
    )

    # Find overlaps: ignore.strand=FALSE ensures strand must match
    overlaps <- findOverlaps(gr_reads, gr_ip, type = "any", ignore.strand = FALSE)

    # Extract overlapping reads and annotations
    merged_reads <- gr_reads[queryHits(overlaps)]
    merged_ip <- gr_ip[subjectHits(overlaps)]
    merged = data.frame(
        seqnames = seqnames(merged_reads),
        read_start = start(merged_reads),
        read_end = end(merged_reads),
        strand = strand(merged_reads),
        anno_start = merged_reads$anno_start,
        anno_end = merged_reads$anno_end,
        startR = merged_reads$startR,
        endR = merged_reads$endR,
        ip_start = start(merged_ip),
        ip_end = end(merged_ip),
        gene_name = merged_reads$gene_name,
        transcript_name = merged_reads$transcript_name,
        fidr = merged_reads$fidr
    ) %>% 
    data.table()

    #Build IP table based on mapping
    ipdb = merged %>%
        filter(ip_start > anno_start & ip_end < anno_end) %>%
        #calculation of relative coordinates of ip positions
        mutate(
            start_ipR = case_when(
                strand == "-" ~ endR - (anno_end - ip_start),
                strand == "+" ~ endR - (ip_end - anno_start)
            ),
            end_ipR = case_when(
                strand == "-" ~ endR - (anno_end - ip_end),
                strand == "+" ~ endR - (ip_start - anno_start)
            )
        ) %>%
        #identification of all fragment associated to ip positions with start_ipR >= args$threshold_dist
        filter(start_ipR >= args$threshold_dist) %>%
        group_by(transcript_name) %>%
        filter(start_ipR > min(start_ipR)) %>%
        ungroup()
    
    ipdb = merged %>%
        filter(ip_start > anno_start & ip_end < anno_end) %>%
        #calculation of relative coordinates of ip positions
        mutate(
            start_ipR = case_when(
            strand == "-" ~ endR - (anno_end - ip_start),
            strand == "+" ~ endR - (ip_end - anno_start)
            ),
            end_ipR = case_when(
            strand == "-" ~ endR - (anno_end - ip_end),
            strand == "+" ~ endR - (ip_start - anno_start)
            )
        ) %>%
        #identification of all fragment associated to ip positions with start_ipR >= 60
        filter(start_ipR >= 60) %>%
        group_by(transcript_name) %>%
        filter(start_ipR > min(start_ipR)) %>%
        ungroup()
    
    #writing retained ip position
    fwrite(
        ipdb %>%
            distinct(seqnames, ip_start, ip_end, strand) %>%
            data.table(),
        file = paste0(ipdb$seqnames[1],".ipdb"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        nThread = as.numeric(args$threads)
    )

    #discard fragments associated to internal priming positions
    filt1_frags = ipdb %>%
        pull(fidr) %>%
        unique()
    reads = reads %>%
        filter(!fidr %in% filt1_frags) %>%
        select(!fidr) %>%
        data.table()
}

#write reads
fst::write_fst(
  reads,
  path = args$output_path,
  compress = 50
)
