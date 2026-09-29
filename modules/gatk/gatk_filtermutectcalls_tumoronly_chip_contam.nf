process GATK_FILTERMUECTCALLS {
    tag "$sampleID"

    cpus 1
    memory  15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_mutect2_somatic.filtered.vcf.gz", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_mutect2_somatic.filtered.vcf.gz.tbi", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.filteringStats.tsv", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi), path(stats), path(contam_table), path(segments), path(read_orientation_model), val(apply_orientation)

    output:
    tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz"), file("*_mutect2_somatic.filtered.vcf.gz.tbi"), emit: mutect2_vcf_tbi
    tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

    script:
    // chip_somatic counterpart of gatk_filtermutectcalls_wes.nf for backgrounds
    // that run contamination estimation (human). Unlike the WES module, the
    // orientation-bias model is not gated on params.ffpe (see finding 6), and
    // --f-score-beta is passed from chip_somatic.config.
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    // Orientation-bias priors are applied only when the reads are NOT duplex
    // consensus (apply_orientation, set by CHIP_SOMATIC_CALLING's caller).
    // Duplex consensus reads all come out F1R2: fgbio 4.0.0 GroupReadsByUmi
    // (Paired) always labels the F1R2 strand /A, and DuplexConsensusCaller
    // builds duplex R1 from the AB R1s. Verified 2026-09-25 on synthetic
    // combined_sample_v2: raw reads ~52/48 F2R1/F1R2, aligned duplex
    // consensus 992,281 F1R2 vs 13 F2R1, so the model flagged 11/15 calls
    // `orientation`. With no F2R1 reads the model can't separate artifacts
    // from real variants in duplex mode.
    orientation_priors = apply_orientation ?
        "--orientation-bias-artifact-priors ${read_orientation_model}" : ''
    f_score_beta = params.f_score_beta ? "--f-score-beta ${params.f_score_beta}" : ''

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" FilterMutectCalls \
    -R ${params.ref_fa} \
    -V ${vcf} \
    --stats ${stats} \
    --contamination-table ${contam_table} \
    --tumor-segmentation ${segments} \
    ${orientation_priors} \
    ${f_score_beta} \
    -O ${sampleID}_mutect2_somatic.filtered.vcf.gz
    """
}
