process FILTER_BEDPE {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID}/merged_sv", pattern: "*.bedpe", mode: 'copy'
    
    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    input:
        tuple val(sampleID), file(sv_genes_cnv_bedpe)
        val(suppl_switch)
    
    output:
        tuple val(sampleID), file("*sv_annotated_high_confidence*.bedpe")

    script:
        if(suppl_switch == "main")
        """
        Rscript ${projectDir}/bin/wgs/filter-bedpe.r \
            --bedpe=${sv_genes_cnv_bedpe} \
            --outfile_highconf=${sampleID}_sv_annotated_high_confidence.bedpe
        """

        else if (suppl_switch == "supplemental")
        """
        Rscript ${projectDir}/bin/wgs/filter-bedpe.r \
            --bedpe=${sv_genes_cnv_bedpe} \
            --outfile_highconf=${sampleID}_sv_annotated_high_confidence_supplemental.bedpe
        """
}
