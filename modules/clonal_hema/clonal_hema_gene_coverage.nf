// clonal_hema per-gene coverage and on-target summary (2026-11-10).
// Inputs come from PICARD_COLLECTHSMETRICS (consensus BAM) and PICARD_COLLECTHSMETRICS_RAW (raw aligned BAM,
// before UMI deduplication), modules/picard/picard_collecthsmetrics_chip.nf: PER_BASE_COVERAGE + HsMetrics.
// Script: modules/clonal_hema/bin/gene_coverage.py, passed as a `path` input so -resume hashes its content
// (same reason as CHIP_GENE_REPORT). Reuses the pipeline's MultiQC container only for Python 3 (stdlib only).

process CLONAL_HEMA_READ_ON_TARGET {
    // Per-read on-target rate (2026-11-11): primary mapped reads (excluding unmapped 0x4, secondary 0x100,
    // QC-fail 0x200, supplementary 0x800) that overlap >= 1 base of a --target_picard interval, divided by all
    // primary mapped reads. Counted on the raw aligned BAM (before UMI deduplication) and the consensus BAM.
    // No index is staged, so samtools view -L filters while streaming and counts each read once (also with
    // overlapping intervals). Also counted with the intervals padded by params.read_on_target_padding bp per side.
    // reads_primary = all primary reads incl. unmapped (-F 0x900): read pairs = reads_primary / 2.
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    input:
    tuple val(sampleID), path(raw_bam), path(cons_bam)
    path(targets)

    output:
    tuple val(sampleID), path("${sampleID}.read_on_target.tsv"), emit: tsv

    script:
    """
    grep -v '^@' ${targets} | awk 'BEGIN{OFS="\\t"}{print \$1,\$2-1,\$3}' > targets.bed
    hdr="sample"; row="${sampleID}"
    for pad in 0 \$(echo "${params.read_on_target_padding ?: ''}" | tr ',' ' '); do
        sfx=\$( [ "\$pad" = 0 ] && echo "" || echo "_pad\${pad}" )
        awk -v p=\$pad 'BEGIN{OFS="\\t"}{s=\$2-p; if (s<0) s=0; print \$1,s,\$3+p}' targets.bed > t\${pad}.bed
        for x in raw dedup; do
            bam=\$( [ "\$x" = raw ] && echo ${raw_bam} || echo ${cons_bam} )
            if [ "\$pad" = 0 ]; then
                hdr="\$hdr\\treads_primary_\$x\\treads_mapped_\$x"
                row="\$row\\t\$(samtools view -c -F 0x900 \$bam)\\t\$(samtools view -c -F 0xB04 \$bam)"
            fi
            hdr="\$hdr\\treads_on_target\${sfx}_\$x"; row="\$row\\t\$(samtools view -c -F 0xB04 -L t\${pad}.bed \$bam)"
        done
    done
    printf "\$hdr\\n\$row\\n" > ${sampleID}.read_on_target.tsv
    """
}

process CLONAL_HEMA_SEX_CHECK {
    // Inferred sex (2026-11-11): primary reads with MAPQ >= 30 (samtools view -c -q 30 -F 0xF04) on the autosomes,
    // chrX and chrY of the raw aligned BAM (before UMI deduplication; indexed, so region queries). Mostly off-target
    // reads for chrY. Contig names are taken from the BAM header (X/Y or chrX/chrY; autosomes = numeric names).
    // The call itself is made in gene_coverage.py from these counts.
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    input:
    tuple val(sampleID), path(raw_bam), path(raw_bai)

    output:
    tuple val(sampleID), path("${sampleID}.sex_check.tsv"), emit: tsv

    script:
    """
    samtools idxstats ${raw_bam} | cut -f1 > contigs.txt
    x=\$(grep -E -x '(chr)?X' contigs.txt | head -1 || true)
    y=\$(grep -E -x '(chr)?Y' contigs.txt | head -1 || true)
    auto=\$(grep -E -x '(chr)?[0-9]+' contigs.txt | tr '\\n' ' ')
    n_auto=\$(samtools view -c -q 30 -F 0xF04 ${raw_bam} \$auto)
    n_x=\$( [ -n "\$x" ] && samtools view -c -q 30 -F 0xF04 ${raw_bam} \$x || echo 0 )
    n_y=\$( [ -n "\$y" ] && samtools view -c -q 30 -F 0xF04 ${raw_bam} \$y || echo 0 )
    printf 'sample\\tsex_autosome_reads\\tsex_chrX_reads\\tsex_chrY_reads\\n%s\\t%s\\t%s\\t%s\\n' ${sampleID} \$n_auto \$n_x \$n_y > ${sampleID}.sex_check.tsv
    """
}

process CLONAL_HEMA_GENE_COVERAGE {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/multiqc:v1.34.dev0_custom'

    publishDir path: { "${params.pubdir}/${sampleID}/stats" }, pattern: "*.tsv", mode: 'copy'

    input:
    tuple val(sampleID), path(dedup_per_base), path(dedup_hsmetrics), path(raw_per_base), path(raw_hsmetrics), path(read_on_target), path(sex_check)
    path(script)

    output:
    tuple val(sampleID), path("${sampleID}.gene_coverage.tsv"), emit: gene_tsv
    tuple val(sampleID), path("${sampleID}.on_target.tsv"), emit: on_target_tsv

    script:
    """
    python3 ${script} sample --sample ${sampleID} \
        --dedup-per-base ${dedup_per_base} --raw-per-base ${raw_per_base} \
        --dedup-hsmetrics ${dedup_hsmetrics} --raw-hsmetrics ${raw_hsmetrics} \
        --read-on-target ${read_on_target} --sex-check ${sex_check} \
        --thresholds ${params.gene_coverage_thresholds}
    """
}

process CLONAL_HEMA_COVERAGE_SUMMARY {

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/multiqc:v1.34.dev0_custom'

    publishDir path: { "${params.pubdir}/coverage" }, pattern: "*.tsv", mode: 'copy'

    input:
    path(gene_tsvs)
    path(on_target_tsvs)
    path(script)

    output:
    path("*_mqc.tsv"), emit: mqc
    path("*.tsv"), emit: tables
    path("gene_coverage_all_samples.tsv"), emit: gene_all
    path("on_target_all_samples.tsv"), emit: on_target_all

    script:
    """
    python3 ${script} summary --gene-tsv ${gene_tsvs} --on-target-tsv ${on_target_tsvs}
    """
}
