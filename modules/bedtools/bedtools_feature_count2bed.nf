process FEATURE_COUNT2BED {
    tag "$sampleID"

    cpus 1
    memory 4.GB     
    time '04:00:00'  
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_peaks_countMatrix.mm10.bed", mode: 'copy'
    
    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(sampleID), path(peak_cnt_matrx)

    output:
    tuple val(sampleID), path("*_peaks_countMatrix.mm10.bed")

    script:
    """
    bash "${moduleDir}/bin/bedtools_feature_count2bed.sh" \
        "${peak_cnt_matrx}" \
        "${sampleID}"
    """
}