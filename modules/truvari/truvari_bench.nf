process TRUVARI_BENCH {
    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "12:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container '/projects/compsci/vmp/USERS/widmas/containers/bcftools_1.22_truvari_5.3.0.img'

    input:
        tuple val(sampleID), file(vcf_tuple)
    
    output:
        tuple val(sampleID), file("${sampleID}_ensemble.vcf.gz"), file("${sampleID}_ensemble.vcf.gz.tbi"), emit: vcf_ensemble
    
    script:

    """
    # Compress VCF files
    bcftools view -Oz -o ${sampleID}_pbsv.vcf.gz ${vcf_tuple[0]}
    bcftools view -Oz -o ${sampleID}_sniffles.vcf.gz ${vcf_tuple[1]}
    
    # Index the compressed VCF files
    bcftools index ${sampleID}_pbsv.vcf.gz
    bcftools index ${sampleID}_sniffles.vcf.gz
    
    # Run truvari comparison
    truvari bench -b ${sampleID}_pbsv.vcf.gz -c ${sampleID}_sniffles.vcf.gz -o truvari_bench \
        -P 0.75 \
        -O 0.5 \
        -p 0.7 \
        -s 20 \
        -S 10000000
    
    # Move the true positive calls to the final output
    mv truvari_bench/tp-base.vcf.gz ${sampleID}_ensemble.vcf.gz
    mv truvari_bench/tp-base.vcf.gz.tbi ${sampleID}_ensemble.vcf.gz.tbi
    """

}