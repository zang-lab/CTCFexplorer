# CTCFexplorer Nextflow Pipeline Example

This pipeline takes a plain `SRR_list.txt`, downloads each accession from SRA, runs basic QC, aligns reads to the genome, cleans BAMs, calls MACS2 peaks, builds a union CTCF site set, and characterizes union sites into high-confidence and constitutive subsets.

## Input

Create a text file with one SRR accession per line:

```text
SRR17619455
SRR17619456
SRR17619457
```

Blank lines and lines starting with `#` are ignored.

## What the pipeline does

1. Downloads each SRR with `prefetch` and `fasterq-dump`.
2. Runs `fastqc` on raw FASTQs.
3. Aligns reads to the supplied Bowtie2 genome index.
4. Cleans BAMs by removing unmapped reads, duplicates, low-MAPQ reads, nonstandard chromosomes, and blacklist overlaps.
5. Calls peaks with MACS2.
6. Builds a union set from summit-centered windows on all peak files.
7. Computes occupancy per union site across all samples.
8. Counts reads on the union sites and converts them to RPKM.
9. Produces:
   - `high_confidence.bed`
   - `constitutive.bed`
   - filtered occupancy matrices
   - filtered RPKM matrices
   - quantile-normalized high-confidence RPKM matrix

## Required software

The workflow assumes these commands are available on the execution nodes:

- `prefetch`
- `fasterq-dump`
- `fastqc`
- `bowtie2`
- `samtools`
- `bedtools`
- `picard`
- `macs3` (or `macs2` if you set `--macs_cmd macs2`)
- `python3`

Python helper scripts also require:

- `pandas`
- `numpy`
- `pysam`

If you prefer modules on UVA HPC, load them before launching Nextflow, or adapt the processes to include your module commands.

## Required references

You must provide:

- `--bowtie2_index`
- `--blacklist_bed`

Example for hg38:

```bash
export NXF_HOME="$(pwd)/.nextflow"

nextflow run main.nf \
  -profile slurm \
  --srr_list /path/to/SRR_list.txt \
  --bowtie2_index /path/to/bowtie2/index/prefix \
  --blacklist_bed /path/to/blacklist.bed \
  --outdir /path/to/results \
  --queue standard \
  --account <slurm_account>
```

Set `NXF_HOME` to any writable directory if your environment does not allow Nextflow to write to the default home location.

## Stub mode

`-stub-run` now writes larger fake FASTQs so downstream tools see nontrivial inputs. By default each stub FASTQ contains `1,000,000` reads. You can override that with:

```bash
--stub_read_count 1000000
```

## Synthetic real-run mode

If you want to test the actual downstream tools without downloading from SRA, the workflow supports a synthetic mode in normal execution. In that mode, `DOWNLOAD_FASTQ` generates real gzipped FASTQs and the rest of the pipeline runs normally.

Required extra parameter:

```bash
--synthetic_mode true \
--synthetic_reference_fasta /path/to/synthetic.fa
```

Optional controls:

```bash
--synthetic_read_count 1000000 \
--synthetic_read_length 75
```

You can create a simple synthetic reference with:

```bash
python3 bin/create_synthetic_reference.py \
  --output synthetic.fa \
  --length 20000
```

## Important defaults

- MACS2: `--nomodel --extsize 146 -g hs -q 0.01`
- Union site width: `150 bp`
- Peak fold-enrichment cutoff for union: `4`
- High-confidence threshold: occupancy in at least `2` samples
- Constitutive threshold: occupancy in `100%` of samples

You can override these:

```bash
--union_width 150 \
--union_fe_cutoff 4 \
--high_confidence_min_samples 2 \
--constitutive_min_fraction 1.0
```

If you specifically want `macs2`, override:

```bash
--macs_cmd macs2
```

## Output layout

Results are written under `--outdir`:

- `01_fastq`
- `02_fastqc`
- `03_align`
- `04_clean_bam`
- `05_peaks`
- `06_union`
- `07_occupancy`
- `08_signal`
- `09_final`
- `pipeline_info`

## Notes

- Sample IDs in this workflow are the SRR accessions themselves.
- The workflow assumes human-style chromosome names (`chr1` ... `chr22`, `chrX`, `chrY`).
- `constitutive.bed` is based on a fraction threshold of total sample count rather than a hard-coded count.
- `nextflow -resume` is the intended way to rerun after failures or parameter changes.
