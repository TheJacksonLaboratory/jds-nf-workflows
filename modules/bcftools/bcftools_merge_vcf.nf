process BCFTOOLS_MERGE_VCF {
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/results/gold_truth_vcf", pattern: '*ALLchr_golden.vcf*', mode:'copy'

    input:
    path(vcf)
    path(tbi)
    val(name)   


    output:
    path("*ALLchr_golden.vcf.gz"), emit: merged_vcf, optional: true
    path("*ALLchr_golden.vcf"), emit: merged_vcf_unzip, optional: true

    script:
    prefix = params.gen_org=='mouse' ? "Mus_musculus.GRCm38" : "Homo_sapiens.GRCh38"
    coverage = params.workflow == "generate_wgs_simreads" ? "10x" : "60x"        
    """
    ls *.vcf.gz > vcfout.list

    bcftools merge -m none --file-list vcfout.list -Oz -o ${prefix}_simVar_${coverage}_${name}chr_golden.vcf.gz

    if [[ ${name} == "FINAL_ALL" ]]; then
        gunzip ${prefix}_simVar_${coverage}_${name}chr_golden.vcf.gz
    fi
    """
}
