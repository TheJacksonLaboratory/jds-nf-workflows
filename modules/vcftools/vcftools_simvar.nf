process VCFTOOLS_SIMVAR {
    tag "intersect"

    cpus = 1
    memory = 10.GB
    time = '23:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/vcftools:0.1.17--g581c231'

    publishDir "${params.pubdir}/results/gold_truth_vcf", pattern: '*intersect.recode*', mode:'copy'


    input:
        path(merged_vcf)
        path(merged_vcf_tbi)

    output:
        path("*_intersect*"), emit: int_vcf
       

    script:
        prefix = params.gen_org=='mouse' ? "Mus_musculus.GRCm38" : "Homo_sapiens.GRCh38"
        target_bed = params.gen_org=='mouse' ? params.target_bed : params.target_sort_bed 
        """
        vcftools \
             --gzvcf ${merged_vcf} \
             --out ${prefix}_simVar_60x_ALLchr_golden_WES_intersect \
             --bed ${target_bed} \
             --recode \
             --recode-INFO-all
        """

}
