
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
parser$add_argument('--annotations', type="character", help='path of isoform annotation file (exons, relative coordinates)')
parser$add_argument('--output_path', type='character', help='path of output reads file')
parser$add_argument("--threads", help = "Number of threads")
parser$add_argument('--threshold_dist', type="character", help='distance of ip position from isoform end')
args = parser$parse_args()

THRESHOLD_IP_DIST = as.numeric(args$threshold_dist)

# Load data
# +++++++++
reads = fst::read_fst(path = args$reads, as.data.table = T)
if(file.size(args$ip_reference) == 0){
  #- no internal priming site on this chromosome...
  ip_pos = data.table(seqnames = character(0))
} else {
  ip_pos = fread(
    file = args$ip_reference, 
    col.names = c("seqnames", "start_ip", "end_ip", "V4", "V5", "strand_ip"),
    nThread = as.numeric(args$threads)
  ) %>%
    filter(seqnames %in% unique(reads$seqnames))
}



# Filtering of internal priming reads (original SCALPEL rule)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  for each isoform, the internal priming (IP) site closest to its 3' end (>= THRESHOLD_IP_DIST)
#  is located on its exons; every fragment x isoform pair lying entirely upstream of it is
#  discarded (fragments there may come from oligo-dT priming on the genomic A-rich stretch)
if(nrow(ip_pos) == 0){
  message("No internal priming positions found in the provided annotation file.")

} else {

    annots = fread(args$annotations) %>%
        dplyr::filter(transcript_name %in% unique(reads$transcript_name)) %>%
        dplyr::distinct(seqnames, start, end, strand, startR, endR, transcript_name) %>%
        data.table()

    gr_annots = GRanges(
        seqnames = annots$seqnames,
        ranges = IRanges(start = annots$start, end = annots$end),
        strand = annots$strand
    )
    gr_ip = GRanges(
        seqnames = ip_pos$seqnames,
        ranges = IRanges(start = ip_pos$start_ip, end = ip_pos$end_ip),
        strand = ip_pos$strand_ip
    )

    # IP sites on the isoforms exons (same strand)
    hits = findOverlaps(gr_annots, gr_ip, type = "any", ignore.strand = FALSE)
    ipdb = cbind(annots[queryHits(hits)], ip_pos[subjectHits(hits), .(start_ip, end_ip)]) %>%
        dplyr::filter(start_ip > start & end_ip < end) %>%
        #- relative coordinates of IP sites (distance to the isoform 3' end)
        dplyr::mutate(
            start_ipR = ifelse(strand == "+", endR - (end_ip - start), endR - (end - start_ip))
        ) %>%
        dplyr::filter(start_ipR >= THRESHOLD_IP_DIST) %>%
        #- closest IP site out of threshold distance for each isoform
        dplyr::group_by(transcript_name) %>%
        dplyr::filter(start_ipR == min(start_ipR)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(transcript_name, strand, start_ip, end_ip) %>%
        data.table()
    message("isoforms with an internal priming site: ", uniqueN(ipdb$transcript_name))

    #writing retained ip position
    fwrite(
        ipdb,
        file = paste0(reads$seqnames[1], ".ipdb"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        nThread = as.numeric(args$threads)
    )

    #- fragments x isoforms lying entirely upstream of the closest IP site (both strands)
    frags = reads[, .(fg_start = min(read_start), fg_end = max(read_end)), by = .(tags, transcript_name)]
    frags_todel = merge(frags, ipdb, by = "transcript_name") %>%
        dplyr::filter((strand == "+" & fg_end < start_ip) | (strand == "-" & fg_start > end_ip)) %>%
        dplyr::distinct(tags, transcript_name) %>%
        data.table()
    message("fragment x isoform pairs discarded (internal priming): ", nrow(frags_todel))

    reads = reads[!frags_todel, on = .(tags, transcript_name)]
}

#write reads
fst::write_fst(
  reads,
  path = args$output_path,
  compress = 50
)
