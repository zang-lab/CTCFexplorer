nextflow.enable.dsl=2

def requireParam(name) {
    if (!params[name]) {
        error "Missing required parameter: --${name}"
    }
}

requireParam('srr_list')
requireParam('bowtie2_index')
requireParam('blacklist_bed')

if (params.synthetic_mode && !params.synthetic_reference_fasta) {
    error "Missing required parameter for synthetic mode: --synthetic_reference_fasta"
}

def srrFile = file(params.srr_list)
if (!srrFile.exists()) {
    error "SRR list not found: ${params.srr_list}"
}

def srrIds = srrFile.readLines()
    .collect { it.trim() }
    .findAll { it && !it.startsWith('#') }
    .unique()

if (!srrIds) {
    error "No SRR accessions found in ${params.srr_list}"
}

process DOWNLOAD_FASTQ {
    tag "$srr"
    label 'download'
    publishDir "${params.outdir}/01_fastq", mode: 'copy'

    input:
    val srr

    output:
    tuple val(srr), path("${srr}*.fastq.gz"), emit: reads
    path("${srr}.download.log"), emit: logs

    script:
    """
    set -euo pipefail

    if [[ "${params.synthetic_mode}" == "true" ]]; then
      generate_synthetic_fastq.py \\
        --sample ${srr} \\
        --reference ${params.synthetic_reference_fasta} \\
        --read-count ${params.synthetic_read_count} \\
        --read-length ${params.synthetic_read_length} \\
        --output ${srr}.fastq.gz \\
        > ${srr}.download.log 2>&1
    else
      mkdir -p cache fastq

      if prefetch --output-directory cache ${srr} > ${srr}.download.log 2>&1; then
        fasterq-dump --threads ${task.cpus} --split-files --skip-technical -O fastq cache/${srr}/${srr}.sra >> ${srr}.download.log 2>&1
      else
        fasterq-dump --threads ${task.cpus} --split-files --skip-technical -O fastq ${srr} >> ${srr}.download.log 2>&1
      fi

      ls fastq/*.fastq >/dev/null
      gzip -f fastq/*.fastq
      mv fastq/*.fastq.gz .
    fi
    """

    stub:
    """
    python3 - <<'PY'
import gzip

srr = "${srr}"
read_count = int("${params.stub_read_count}")
sequence = "ACGT" * 25
quality = "F" * len(sequence)

with gzip.open(f"{srr}.fastq.gz", "wt") as out:
    for i in range(1, read_count + 1):
        out.write(f"@{srr}.{i}\\n{sequence}\\n+\\n{quality}\\n")
PY
    : > ${srr}.download.log
    """
}

process FASTQC_RAW {
    tag "$srr"
    label 'qc'
    publishDir "${params.outdir}/02_fastqc", mode: 'copy'

    input:
    tuple val(srr), path(reads)

    output:
    path("*_fastqc.zip"), emit: zips
    path("*_fastqc.html"), emit: html

    script:
    def readArgs = reads.collect { it.getName() }.join(' ')
    """
    set -euo pipefail
    fastqc --threads ${task.cpus} ${readArgs}
    """

    stub:
    """
    for f in ${reads.collect { it.getName() }.join(' ')}; do
      base=\$(basename "\$f" .fastq.gz)
      : > \${base}_fastqc.zip
      : > \${base}_fastqc.html
    done
    """
}

process ALIGN_AND_CLEAN {
    tag "$srr"
    label 'align'
    publishDir "${params.outdir}/03_align", mode: 'copy', pattern: "${srr}.bam"
    publishDir "${params.outdir}/04_clean_bam", mode: 'copy', pattern: "${srr}_clean.bam*"
    publishDir "${params.outdir}/04_clean_bam", mode: 'copy', pattern: "${srr}_clean.flagstat.txt"
    publishDir "${params.outdir}/04_clean_bam", mode: 'copy', pattern: "${srr}_dup_metrics.txt"

    input:
    tuple val(srr), path(reads)

    output:
    tuple val(srr), path("${srr}.bam"), path("${srr}_clean.bam"), path("${srr}_clean.bam.bai"), path("${srr}_clean.flagstat.txt"), emit: cleaned

    script:
    def alignArgs = reads.size() == 2
        ? "-1 ${reads[0].getName()} -2 ${reads[1].getName()}"
        : "-U ${reads.collect { it.getName() }.join(',')}"
    """
    set -euo pipefail

    bowtie2 -x ${params.bowtie2_index} ${alignArgs} -p ${task.cpus} \\
      | samtools view -@ ${task.cpus} -bS -o ${srr}.bam

    samtools view -@ ${task.cpus} -b -F 4 ${srr}.bam \\
      | samtools sort -@ ${task.cpus} -o ${srr}_tmp_sorted.bam

    picard MarkDuplicates \\
      I=${srr}_tmp_sorted.bam \\
      O=${srr}_tmp_marked.bam \\
      M=${srr}_dup_metrics.txt \\
      REMOVE_DUPLICATES=false \\
      VALIDATION_STRINGENCY=SILENT

    samtools view -h ${srr}_tmp_marked.bam \\
      | awk 'substr(\$0,1,1)=="@" || \$3 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)\$/' \\
      | samtools view -@ ${task.cpus} -b -q ${params.mapq} -F 1024 - \\
      | bedtools intersect -abam - -b ${params.blacklist_bed} -v \\
      > ${srr}_clean.bam

    samtools index ${srr}_clean.bam
    samtools flagstat ${srr}_clean.bam > ${srr}_clean.flagstat.txt
    """

    stub:
    """
    : > ${srr}.bam
    : > ${srr}_clean.bam
    : > ${srr}_clean.bam.bai
    printf "0 + 0 in total\\n" > ${srr}_clean.flagstat.txt
    printf "LIBRARY\\tUNPAIRED_READS_EXAMINED\\n" > ${srr}_dup_metrics.txt
    """
}

process CALL_PEAKS {
    tag "$srr"
    label 'peak'
    publishDir "${params.outdir}/05_peaks", mode: 'copy'

    input:
    tuple val(srr), path(raw_bam), path(clean_bam), path(clean_bai), path(flagstat)

    output:
    path("${srr}_peaks.narrowPeak"), emit: narrowpeaks
    path("${srr}_summits.bed"), emit: summits
    path("${srr}_qc.tsv"), emit: peak_qc

    script:
    """
    set -euo pipefail

    ${params.macs_cmd} callpeak \\
      --nomodel \\
      --extsize ${params.macs2_extsize} \\
      -g ${params.macs2_genome_size} \\
      -B \\
      --keep-dup 1 \\
      -q ${params.macs2_qvalue} \\
      --SPMR \\
      -n ${srr} \\
      -t ${clean_bam}

    total_reads=\$(samtools view -c ${raw_bam})
    clean_reads=\$(samtools view -c ${clean_bam})
    frip_reads=\$(bedtools intersect -abam ${clean_bam} -b ${srr}_peaks.narrowPeak -u | samtools view -c)
    frip=\$(awk -v frip=\$frip_reads -v total=\$clean_reads 'BEGIN { if (total > 0) printf "%.4f", frip / total; else printf "0" }')
    peak_count=\$(wc -l < ${srr}_peaks.narrowPeak)

    printf "sample\\ttotal_reads\\tclean_reads\\tfrip\\tpeak_count\\n" > ${srr}_qc.tsv
    printf "${srr}\\t%s\\t%s\\t%s\\t%s\\n" "\$total_reads" "\$clean_reads" "\$frip" "\$peak_count" >> ${srr}_qc.tsv
    """

    stub:
    """
    cat <<'EOF' > ${srr}_peaks.narrowPeak
    chr1	100	250	${srr}_peak1	1000	.	10	20	5	75
    chr1	500	650	${srr}_peak2	800	.	8	16	4	60
    EOF
    cat <<'EOF' > ${srr}_summits.bed
    chr1	175	176
    chr1	560	561
    EOF
    printf "sample\\ttotal_reads\\tclean_reads\\tfrip\\tpeak_count\\n${srr}\\t100\\t80\\t0.2500\\t2\\n" > ${srr}_qc.tsv
    """
}

process BUILD_UNION {
    label 'aggregate'
    publishDir "${params.outdir}/06_union", mode: 'copy'

    input:
    path(narrowpeaks)

    output:
    path("all_peaks.narrowPeak"), emit: all_peaks
    path("union_sites.bed"), emit: union_bed

    script:
    def peakArgs = narrowpeaks.collect { it.getName() }.join(' ')
    """
    set -euo pipefail

    cat ${peakArgs} > all_peaks.narrowPeak

    create_union.py \\
      --inputs ${peakArgs} \\
      --output union_sites.bed \\
      --width ${params.union_width} \\
      --fe-cutoff ${params.union_fe_cutoff}
    """

    stub:
    """
    cat ${narrowpeaks.collect { it.getName() }.join(' ')} > all_peaks.narrowPeak
    cat <<'EOF' > union_sites.bed
    chr1	100	250	union_0000001	0	+
    chr1	500	650	union_0000002	0	+
    EOF
    """
}

process BUILD_OCCUPANCY_MATRIX {
    label 'aggregate'
    publishDir "${params.outdir}/07_occupancy", mode: 'copy'

    input:
    path union_bed
    path narrowpeaks

    output:
    path("union_occupancy_matrix.tsv"), emit: matrix
    path("union_occupancy_scores.tsv"), emit: scores

    script:
    def peakArgs = narrowpeaks.collect { it.getName() }.join(' ')
    """
    set -euo pipefail

    build_occupancy_matrix.py \\
      --union ${union_bed.getName()} \\
      --peaks ${peakArgs} \\
      --output-matrix union_occupancy_matrix.tsv \\
      --output-scores union_occupancy_scores.tsv \\
      --fe-cutoff ${params.union_fe_cutoff} \\
      --sample-width ${params.occupancy_width}
    """

    stub:
    def sampleNames = narrowpeaks.collect { it.getName().replace('_peaks.narrowPeak', '') }
    def header = (['Union_ID'] + sampleNames).join('\\t')
    def ones = sampleNames.collect { '1' }.join('\\t')
    """
    cat <<'EOF' > union_occupancy_matrix.tsv
    ${header}
    union_0000001	${ones}
    union_0000002	${ones}
    EOF
    cat <<'EOF' > union_occupancy_scores.tsv
    chr	start	end	Union_ID	score	strand	occupancy_score
    chr1	100	250	union_0000001	0	+	${sampleNames.size()}
    chr1	500	650	union_0000002	0	+	${sampleNames.size()}
    EOF
    """
}

process COUNT_SIGNAL_ON_UNION {
    tag "$srr"
    label 'signal'
    publishDir "${params.outdir}/08_signal/per_sample", mode: 'copy'

    input:
    path union_bed
    tuple val(srr), path(raw_bam), path(clean_bam), path(clean_bai), path(flagstat)

    output:
    path("${srr}.signal.tsv"), emit: signal_tables

    script:
    """
    set -euo pipefail

    count_signal_on_union.py \\
      --union ${union_bed.getName()} \\
      --bam ${clean_bam.getName()} \\
      --sample ${srr} \\
      --output ${srr}.signal.tsv
    """

    stub:
    """
    cat <<'EOF' > ${srr}.signal.tsv
    Union_ID	${srr}
    union_0000001	10.0
    union_0000002	20.0
    EOF
    """
}

process MERGE_SIGNAL_MATRIX {
    label 'aggregate'
    publishDir "${params.outdir}/08_signal", mode: 'copy'

    input:
    path signal_tables

    output:
    path("union_RPKM.tsv"), emit: matrix

    script:
    def signalArgs = signal_tables.collect { it.getName() }.join(' ')
    """
    set -euo pipefail

    merge_signal_matrix.py \\
      --inputs ${signalArgs} \\
      --output union_RPKM.tsv
    """

    stub:
    """
    python3 ${projectDir}/bin/merge_signal_matrix.py --inputs ${signal_tables.collect { it.getName() }.join(' ')} --output union_RPKM.tsv
    """
}

process CHARACTERIZE_UNION {
    label 'aggregate'
    publishDir "${params.outdir}/09_final", mode: 'copy'

    input:
    path union_bed
    path occupancy_matrix
    path occupancy_scores
    path union_rpkm
    val sample_count

    output:
    path("high_confidence.bed")
    path("constitutive.bed")
    path("union_occupancy_high_confidence.tsv")
    path("union_occupancy_constitutive.tsv")
    path("union_RPKM_high_confidence.tsv")
    path("union_RPKM_constitutive.tsv")
    path("union_RPKM_high_confidence_QN.tsv")

    script:
    """
    set -euo pipefail

    characterize_union.py \\
      --union ${union_bed.getName()} \\
      --occupancy-matrix ${occupancy_matrix.getName()} \\
      --occupancy-scores ${occupancy_scores.getName()} \\
      --signal-matrix ${union_rpkm.getName()} \\
      --sample-count ${sample_count} \\
      --high-confidence-min-samples ${params.high_confidence_min_samples} \\
      --constitutive-min-fraction ${params.constitutive_min_fraction}
    """

    stub:
    """
    python3 ${projectDir}/bin/characterize_union.py \\
      --union ${union_bed.getName()} \\
      --occupancy-matrix ${occupancy_matrix.getName()} \\
      --occupancy-scores ${occupancy_scores.getName()} \\
      --signal-matrix ${union_rpkm.getName()} \\
      --sample-count ${sample_count} \\
      --high-confidence-min-samples ${params.high_confidence_min_samples} \\
      --constitutive-min-fraction ${params.constitutive_min_fraction}
    """
}

workflow {
    Channel.fromList(srrIds).set { srr_ch }
    Channel.value(srrIds.size()).set { sample_count_ch }

    downloaded = DOWNLOAD_FASTQ(srr_ch)

    FASTQC_RAW(downloaded.reads)

    aligned = ALIGN_AND_CLEAN(downloaded.reads)

    peaks = CALL_PEAKS(aligned.cleaned)

    union = BUILD_UNION(peaks.narrowpeaks.collect())

    occupancy = BUILD_OCCUPANCY_MATRIX(union.union_bed, peaks.narrowpeaks.collect())
    signal = COUNT_SIGNAL_ON_UNION(union.union_bed, aligned.cleaned)
    merged = MERGE_SIGNAL_MATRIX(signal.signal_tables.collect())

    CHARACTERIZE_UNION(
        union.union_bed,
        occupancy.matrix,
        occupancy.scores,
        merged.matrix,
        sample_count_ch
    )
}
