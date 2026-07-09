process SAMTOOLS_FASTQ {
    tag "$sampleID"

    cpus 4
    memory  20.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*.fastq.gz", mode: 'copy'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.fastq.gz"), emit: fastq

    script:
    // samtools fastq is much more memory-efficient than bedtools and streams directly to gzip
    // note, for PE reads: '-0 /dev/null' discards unpaired reads, and '-s /dev/null' discards singletons
    // Change the `-n` to `-N` flag to append /1 and /2 suffixes if needed.
    if (params.read_type == 'PE') {
        """
        samtools fastq \
        -@ ${task.cpus - 1} \
        -1 >(gzip -c > ${sampleID}.R1.fastq.gz) \
        -2 >(gzip -c > ${sampleID}.R2.fastq.gz) \
        -0 /dev/null \
        -s /dev/null \
        -n \
        ${bam}
        """
    } else {
        """
        samtools fastq \
        -@ ${task.cpus - 1} \
        -n \
        ${bam} | gzip -c > ${sampleID}.R1.fastq.gz
        """
    }
}
