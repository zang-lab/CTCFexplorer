#!/bin/bash

module load gcc
module load macs2
module load samtools
module load bedtools

# Single-end alignment
bowtie2 -x $INDEX -U SRR_id.fastq -p 8 | \
  samtools view -Sb -o $GSM_id.bam

# Paired-end alignment
bowtie2 -x $INDEX -1 $SRR_id_1.fastq -2 $SRR_id_2.fastq -p 8 | \
  samtools view -Sb -o $GSM_id.bam

# Define input/output variables within the script
BAM_IN="GSM_id.bam"
CLEAN_OUT="GSM_id_clean.bam"
TMP_SORTED="GSM_id_tmp_sorted.bam"
TMP_MARKED="GSM_id_tmp_marked_dup.bam"
DUP_METRICS="GSM_id_dup_metrics.txt"
BLACKLIST="hg38-blacklist.bed"

# Filter out unmapped reads (-F 4) and sort
samtools view -@ 8 -b -F 4 "$BAM_IN" \
  | samtools sort -@ 8 -o "$TMP_SORTED"

# Mark duplicates with Picard
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I="$TMP_SORTED" \
  O="$TMP_MARKED" \
  M="$DUP_METRICS" \
  REMOVE_DUPLICATES=false

# Remove the now-unused sorted file
rm "$TMP_SORTED"

# Filter to keep only standard autosomes (chr1–22), chrX, chrY:
#    - Keep header lines
#    - Keep alignment lines only if the 3rd field is chr1–22 or chrX or chrY
#    Then remove duplicates (flag 1024), remove MAPQ < 30,
#    and exclude blacklist regions.
samtools view -h "$TMP_MARKED" \
  | awk 'substr($0,1,1)=="@" || $3 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' \
  | samtools view -b -q 30 -F 1024 - \
  | bedtools intersect -abam - -b "$BLACKLIST" -v \
  > "$CLEAN_OUT"

# Remove the temporary marked-dup file
rm "$TMP_MARKED"

# Index the final clean BAM
samtools index "$CLEAN_OUT"

# Call peaks
macs2 callpeak --nomodel --extsize 146 -g hs -B --keep-dup 1 -q 0.01 --SPMR -n GSM_id -t GSM_id_clean.bam