# CTCFexplorer Pipeline

This repository contains the core steps of the ChIP-seq processing and analysis workflow used in **CTCFexplorer**.


Step Overview

`s1_peak_calling`
Align ChIP-seq reads to the genome using **Bowtie2** and call peaks using **MACS2**. This step generates peak files and signal tracks for downstream analysis.


`s2_summit_interval`
- Collect all peaks from all samples into a unified all_peaks.bed file.
- Retain only peaks with fold enrichment > 4.
- Plot the distribution of summit-to-summit distances to determine a cutoff.
- This cutoff distinguishes between the same binding site vs. distinct binding sites.


`s3_create_union`
Merge overlapping summit-centered extended peak regions using the distance cutoff to define a set of **union binding sites** shared across samples.


`s4_occupancy_matrix`
Generate a **binary occupancy matrix**:
- **Rows** = union binding sites  
- **Columns** = samples  
- **Values** = `1` if a sample contains a peak in that site, `0` otherwise


`s5_count_reads_on_union`
- Count reads overlapping each union site in each sample.
- Normalize counts to RPKM values.
- Output a **RPKM matrix** analogous to the occupancy matrix  
  (same rows and columns, but with RPKM values instead of binary indicators).
