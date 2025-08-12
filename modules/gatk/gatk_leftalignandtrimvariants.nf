process GATK_LEFTALIGNANDTRIMVARIANTS {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*_splitAndPassOnly.vcf*", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(tbi)
    val(type)

    output:
    tuple val(sampleID), file("*.split.vcf"), emit: split_vcf
    tuple val(sampleID), file("*_splitAndPassOnly.vcf"), emit: vcf, optional: true
    
    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" LeftAlignAndTrimVariants \
        -R ${params.mt_fasta} \
        -V ${vcf} \
        -O ${sampleID}.split.vcf \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac

    if [[ "${type}" == "pass-only" ]]; then
        gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" SelectVariants \
            -V ${sampleID}.split.vcf \
            -O ${sampleID}_splitAndPassOnly.vcf \
            --exclude-filtered
    fi

    """
}
