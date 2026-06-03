process VCF_SORT {
    tag "sort_merged_vcf"

    cpus = 1
    memory = 10.GB
    time = '23:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/vcftools:0.1.17--g581c231'

    publishDir "${params.pubdir}/results/gold_truth_vcf", pattern: '*sorted.vcf.gz*', mode:'copy'

    input:
        path(merged_vcf)

    output:
        path("*sorted.vcf.gz"), emit: vcf
        path("*sorted.vcf.gz.tbi"), emit: tbi


    script:
        prefix="\$(echo -e ${merged_vcf} | sed 's/.vcf.gz//g')"  
        """
        vcf-sort ${merged_vcf} >${prefix}_sorted.vcf
        bgzip ${prefix}_sorted.vcf
        tabix -f -p vcf ${prefix}_sorted.vcf.gz
        """

} 
            
