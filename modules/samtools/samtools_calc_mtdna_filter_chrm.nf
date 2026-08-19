process CALC_MTDNA_FILTER_CHRM {
    tag "$sampleID"

    cpus 4
    memory 20.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*_mtDNA_Content.txt", mode: 'copy'

    input:
    tuple val(sampleID), path(rmdup_bam_file), path(rmdup_bai_file)

    output:
    tuple val(sampleID), path("*.sorted.rmDup.rmChrM.bam"), emit: rmChrM_bam
    tuple val(sampleID), path("*.sorted.rmDup.rmChrM.bam.bai"), emit: rmChrM_bai
    tuple val(sampleID), path("*_mtDNA_Content.txt"), emit: mtdna_log

    script:
    // Get Mitochondrial and total read counts, calculate %mtDNA and filter Mitochondrial Reads from bam file

    mt_name = params.gen_org == 'mouse' ?  'MT' : 'chrM'

    """
    bash "${moduleDir}/bin/samtools_calc_mtdna_filter_chrm.sh" \
        "${rmdup_bam_file}" \
        "${sampleID}" \
        "${mt_name}" \
        "${task.cpus}"
    """
}
