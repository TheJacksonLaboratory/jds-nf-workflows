process GATK_MUTECT2_MT {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam), path(bai)
    val(type) // mt or shifted_mt
    val(interval)

    output:
    tuple val(sampleID), path("*.vcf.gz"), emit: vcf
    tuple val(sampleID), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(sampleID), path("*.stats"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" Mutect2 \
    -R ${params.mt_fasta} \
    -I ${bam} \
    --read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --read-filter MateUnmappedAndUnmappedReadFilter \
    -O ${sampleID}_${type}.vcf.gz \
    -L ${interval} \
    --annotation StrandBiasBySample \
    --mitochondria-mode \
    --max-reads-per-alignment-start 75 \
    --max-mnp-distance 0
    """
}
