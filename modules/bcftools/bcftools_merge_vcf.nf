process BCFTOOLS_MERGE_VCF {
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    path(vcf)
    path(tbi)

    output:
    path("*ALLchr_golden.vcf.gz"), emit: merged_vcf

    script:
     
    """
    ls *.vcf.gz > vcfout.list

    bcftools merge -m none --file-list vcfout.list -Oz -o ${params.sampleID}_${params.gen_org}_simVar_${params.coverage}x_ALLchr_golden.vcf.gz
    """
}
