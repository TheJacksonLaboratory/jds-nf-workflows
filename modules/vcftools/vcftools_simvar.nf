process VCFTOOLS_SIMVAR {
    tag "intersect"

    cpus 1
    memory 10.GB
    time '23:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/vcftools:0.1.17--g581c231'

    publishDir path: { "${params.pubdir}/gold_truth_vcf" }, pattern: '*recode*.vcf.gz', mode:'copy'

    input:
        path(merged_vcf)
        path(merged_vcf_tbi)
        path(target_bed)

    output:
        path("*recode*.vcf.gz"), emit: int_vcf

    script:
        """
        vcftools \
             --gzvcf ${merged_vcf} \
             --out ${params.sampleID}_simVar_${params.coverage}x_ALLchr_golden_WEStargetIntersect \
             --bed ${target_bed} \
             --recode \
             --recode-INFO-all
        bgzip ${params.sampleID}_simVar_${params.coverage}x_ALLchr_golden_WEStargetIntersect.recode.vcf
        """
}
