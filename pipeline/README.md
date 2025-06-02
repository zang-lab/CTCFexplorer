CTCFexplorer Pipeline

This repository contains the core steps of the ChIP-seq processing and analysis workflow used in CTCFexplorer.

⸻

s1_peak_calling

Align ChIP-seq reads to the genome using Bowtie2 and call peaks using MACS2. This step generates peak files and signal tracks for downstream analysis.

⸻

s2_summit_interval

Aggregate peaks from all samples into a single all_peaks.bed file. Retain only peaks with fold enrichment > 4. Plot the distribution of summit-to-summit distances to define a cutoff that distinguishes between peaks representing the same binding site and distinct binding sites.

⸻

s3_create_union

Merge overlapping summit-centered extended peak regions (based on the cutoff) to create a set of union binding sites.

⸻

s4_occupancy_matrix

Build a binary occupancy matrix: rows represent union binding sites, columns represent samples, values indicate whether a peak is present (1) or absent (0) in each sample at each union site.

⸻

s5_count_reads_on_union

Quantify ChIP-seq read counts at each union site for all samples. Convert counts to RPKM values and merge results into a quantitative matrix, analogous in structure to the occupancy matrix.