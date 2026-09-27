process BCFTOOLS_STATS {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.bcftools_stats.txt", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.filter_breakdown.txt", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi)

    output:
    tuple val(sampleID), path("*.bcftools_stats.txt"), emit: stats
    tuple val(sampleID), path("*.filter_breakdown.txt"), emit: filter_breakdown

    script:
    """
    bcftools stats -s - ${vcf} > ${sampleID}.bcftools_stats.txt

    echo -e "count\\tFILTER" > ${sampleID}.filter_breakdown.txt
    bcftools query -f '%FILTER\\n' ${vcf} | sort | uniq -c | awk '{print \$1"\\t"\$2}' >> ${sampleID}.filter_breakdown.txt
    """
}
