process MUTSERVE {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/genepi/mtdna-server-2:v2.1.16'

    publishDir "${params.pubdir}/${sampleID + '/mt_callers'}", pattern: "*.mutserve.norm.vcf.gz*", mode:'copy'

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("*.mutserve.norm.vcf.gz"), path("*.mutserve.norm.vcf.gz.tbi"), emit: vcf_tbi
    
    script:
    
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    java -Xmx${my_mem}G -jar /opt/mutserve/mutserve.jar \
        call \
        --level ${params.detection_limit} \
        --reference ${params.mt_fasta } \
        --mapQ ${params.mapQ} \
        --baseQ ${params.baseQ} \
        --output ${sampleID}.mutserve.vcf.gz \
        --contig-name ${params.mt_contig_name} \
        --no-ansi \
        --strand-bias 1.6 \
        --write-raw \
        ${bam} 

    bcftools norm \
        -m-any \
        -f ${params.mt_fasta} \
        -o ${sampleID}.mutserve.norm.vcf.gz -Oz \
        ${sampleID}.mutserve.vcf.gz

    tabix -f ${sampleID}.mutserve.norm.vcf.gz
    """
}


