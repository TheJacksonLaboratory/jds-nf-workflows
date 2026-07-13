process FINAL_CALC_FRIP {
    tag "$sampleID"

    cpus  1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*_Fraction_reads_in_peak.txt", mode: 'copy'
    
    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), file(processed_bams), file(reads_peaks_bams)

    output:
    tuple val(sampleID), file("*_Fraction_reads_in_peak.txt")

    script:
    // Calculate fraction of reads in peak
    """
    bash "${moduleDir}/bin/samtools_final_calc_frip.sh" \
        "${processed_bams[0]}" \
        "${reads_peaks_bams[0]}" \
        "${sampleID}"
    """
}
