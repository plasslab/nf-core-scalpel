#!/usr/bin/env bash

# Generation of the light test dataset (-profile test)
# ====================================================
# Subset of subsample_data/ (GSE104556, 2 samples) restricted to a gene-dense
# window of chr17, such that the whole pipeline runs in a few minutes locally.
#
# Usage: bash tests/data/make_test_data.sh [SRC_DIR] [OUT_DIR]
# Requirements: samtools, seqkit, python3 (scipy)

set -euo pipefail

SRC=${1:-subsample_data}
OUT=${2:-tests/data/light}
CHR=chr17
START=14000000
END=16000000
SAMPLES="SRR6129050 SRR6129051"

mkdir -p ${OUT}/annotations ${OUT}/reads ${OUT}/barcodes
TMP=$(mktemp -d)
trap "rm -rf ${TMP}" EXIT


#- annotation: genes fully included in the window...
echo "GTF subsetting..."
gzip -dc ${SRC}/annotations/reduced_GENCODE.vM21.annotation.gtf.gz | \
    awk -F'\t' -v c=${CHR} -v s=${START} -v e=${END} '
        $1==c && $3=="gene" && $4>=s && $5<=e { match($9, /gene_id "[^"]+"/); ids[substr($9, RSTART, RLENGTH)] = 1 }
        { match($9, /gene_id "[^"]+"/); if (substr($9, RSTART, RLENGTH) in ids) print }
    ' | gzip -n > ${OUT}/annotations/genes.gtf.gz

#- transcriptome: transcripts of the selected genes...
echo "Transcriptome subsetting..."
gzip -dc ${OUT}/annotations/genes.gtf.gz | \
    awk -F'\t' '$3=="transcript" { match($9, /transcript_id "[^"]+"/); print substr($9, RSTART+15, RLENGTH-16) }' | \
    sort -u > ${TMP}/transcripts.txt
seqkit grep --id-regexp '^([^|]+)\|' -f ${TMP}/transcripts.txt \
    ${SRC}/annotations/gencode.vM21.transcripts.fa.gz | gzip -n > ${OUT}/annotations/transcripts.fa.gz

#- genome (salmon decoys only): window sequence, keeping the chromosome name...
echo "Genome subsetting..."
seqkit grep -p ${CHR} ${SRC}/annotations/GRCm38.subset.genome.fa.gz | \
    seqkit subseq -r $((START+1)):${END} | \
    seqkit replace -p '.*' -r ${CHR} | gzip -n > ${OUT}/annotations/genome.fa.gz

#- internal priming reference: sites in the window...
echo "Internal priming sites subsetting..."
gzip -dc ${SRC}/annotations/mm10_polya.track.gz | \
    awk -F'\t' -v c=${CHR} -v s=${START} -v e=${END} '$1==c && $2>=s && $3<=e' | \
    gzip -n > ${OUT}/annotations/polya.track.gz

for SAMPLE in ${SAMPLES}; do
    echo "Sample ${SAMPLE}..."
    OUTS=${OUT}/cellranger/${SAMPLE}/outs
    mkdir -p ${OUTS}/filtered_feature_bc_matrix

    #- BAM: reads in the window...
    samtools view -b -o ${OUTS}/possorted_genome_bam.bam \
        ${SRC}/cellranger_samples/${SAMPLE}/outs/possorted_genome_bam.bam ${CHR}:${START}-${END}
    samtools index ${OUTS}/possorted_genome_bam.bam

    #- FASTQ: read pairs of the BAM reads (salmon pseudobulk)...
    samtools view ${OUTS}/possorted_genome_bam.bam | cut -f1 | sort -u > ${TMP}/${SAMPLE}.qnames
    for R in R1 R2; do
        seqkit grep -f ${TMP}/${SAMPLE}.qnames \
            ${SRC}/reads/${SAMPLE}_1000000_S1_L001_${R}_001.fastq.gz | \
            gzip -n > ${OUT}/reads/${SAMPLE}_${R}.fastq.gz
    done

    #- DGE matrix: genes of the window, all cells...
    python3 - ${SRC}/cellranger_samples/${SAMPLE}/outs/filtered_feature_bc_matrix ${OUTS}/filtered_feature_bc_matrix ${OUT}/annotations/genes.gtf.gz <<'PY'
import sys, gzip, re, shutil
import scipy.io
src, out, gtf = sys.argv[1:4]
#- gene ids (unversioned) of the window...
ids = set()
with gzip.open(gtf, "rt") as f:
    for line in f:
        m = re.search(r'gene_id "([^".]+)', line)
        if m: ids.add(m.group(1))
with gzip.open(f"{src}/features.tsv.gz", "rt") as f:
    features = [l.rstrip("\n") for l in f]
keep = [i for i, l in enumerate(features) if l.split("\t")[0] in ids]
mtx = scipy.io.mmread(gzip.open(f"{src}/matrix.mtx.gz")).tocsr()[keep, :]
with gzip.open(f"{out}/matrix.mtx.gz", "wb") as f:
    scipy.io.mmwrite(f, mtx.tocoo(), field="integer")
with gzip.open(f"{out}/features.tsv.gz", "wt") as f:
    f.write("\n".join(features[i] for i in keep) + "\n")
shutil.copy(f"{src}/barcodes.tsv.gz", f"{out}/barcodes.tsv.gz")
print(f"  DGE: {len(keep)} genes x {mtx.shape[1]} cells")
PY

    #- barcodes whitelist...
    cp ${SRC}/test_barcodes/${SAMPLE}_barcodes.txt ${OUT}/barcodes/
done

#- samplesheets (paths relative to the samplesheet directory)...
echo "sample,fastq_1,fastq_2,cranger_path" > ${OUT}/samplesheet.csv
echo "sample,barcodes" > ${OUT}/barcodes_whitelist.csv
for SAMPLE in ${SAMPLES}; do
    echo "${SAMPLE},reads/${SAMPLE}_R1.fastq.gz,reads/${SAMPLE}_R2.fastq.gz,cellranger/${SAMPLE}" >> ${OUT}/samplesheet.csv
    echo "${SAMPLE},barcodes/${SAMPLE}_barcodes.txt" >> ${OUT}/barcodes_whitelist.csv
done

echo "Done."
du -sh ${OUT}
