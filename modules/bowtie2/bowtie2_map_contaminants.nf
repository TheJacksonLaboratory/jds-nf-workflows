process BOWTIE_MAP_CONTAMINANTS {
    tag "$sampleID"

    cpus 16
    memory 30.GB
    time '48:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*_bowtie.log", mode: 'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "filtered.*.stats" , mode: 'copy'

    container 'quay.io/biocontainers/bowtie2:2.4.5--py36hfca12d5_2' 

    input:
    tuple val(sampleID), path(reads)
    path index
    val contaminant_type

    output:
    tuple val(sampleID), path("*sam")                               , emit: bam
    tuple val(sampleID), path('*.filter.unmapped.contaminant.fastq'), emit: unmapped
    path "filtered.*.stats"                                     , emit: stats
    tuple val(sampleID), path("*_bowtie.log"), emit: bowtie_log

    script:
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    bowtie2 \\
        --threads ${task.cpus} \\
        --very-sensitive-local \\
        -k 1 \\
        -x \$INDEX \\
        --un ${sampleID}.${contaminant_type}.filter.unmapped.contaminant.fastq \\
        ${reads} \\
        -S ${sampleID}.filter.contaminant.sam > ${sampleID}.${contaminant_type}_bowtie.log 2>&1

    # extracting number of reads from bowtie logs
    awk -v type=${contaminant_type} 'BEGIN{tot=0} {if(NR==4 || NR == 5){tot += \$1}} END {print "\\""type"\\": "tot }' ${sampleID}.${contaminant_type}_bowtie.log | tr -d , > filtered.${sampleID}_${contaminant_type}.stats

    """

}
