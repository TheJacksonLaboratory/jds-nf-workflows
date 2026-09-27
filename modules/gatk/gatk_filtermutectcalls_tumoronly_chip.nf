process GATK_FILTERMUECTCALLS {
    tag "$sampleID"

    cpus 1
    memory  15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_mutect2_somatic.filtered.vcf.gz", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.filteringStats.tsv", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi), path(stats), path(read_orientation_model), val(apply_orientation)

    output:
    tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz"), file("*_mutect2_somatic.filtered.vcf.gz.tbi"), emit: mutect2_vcf_tbi
    tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

    script:
    // Same GATK version/command shape as gatk_filtermutectcalls_wes.nf, minus the
    // contamination-table/tumor-segmentation inputs — for chip_somatic backgrounds
    // (e.g. B6) where contamination estimation is not run at all (see finding 3,
    // germline_resource_design.md). Orientation-bias modeling is not gated on
    // `params.ffpe` for chip_somatic (finding 6), unlike the WES module.
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
    // No-op unless the caller's config defines this param.
    f_score_beta = params.containsKey('f_score_beta') ?
        "--f-score-beta ${params.f_score_beta}" : ''

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" FilterMutectCalls \
    -R ${params.ref_fa} \
    -V ${vcf} \
    --stats ${stats} \
    ${orientation_priors} \
    ${f_score_beta} \
    -O ${sampleID}_mutect2_somatic.filtered.vcf.gz
    """
}
