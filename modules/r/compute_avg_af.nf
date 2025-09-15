process COMPUTE_AVG_AF {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time "4:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    input:
        tuple val(sampleID), path(vcf)

    output:
        tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    
    """
        Rscript ${projectDir}/bin/mitochondria_variant_calling/compute_avgAF.r \
        --vcf ${vcf} \
        --sampleID ${sampleID}_avgAF \
        --output ${sampleID}_mtdna_mergedCallers_avgAF.vcf
    """
}
