process ANNOTATE_GENES_SV {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID}/merged_sv", pattern:"*.MDLS_sv_annotated_genes*.bed", mode:'copy'

    input:
        tuple val(sampleID), file(annot_sv_bedpe)
        val(suppl_switch)

    output:
        tuple val(sampleID), file("*.MDLS_sv_annotated_genes*.bed"), emit: annot_sv_genes_bedpe

    script:
    if (suppl_switch == "main") {
        """
        Rscript ${projectDir}/bin/wgs/annotate-bedpe-with-genes.r \
            --ensembl=${params.ensemblUniqueBed} \
            --bedpe=${annot_sv_bedpe} \
            --out_file=${sampleID}.MDLS_sv_annotated_genes.bed
        """
    } else if (suppl_switch == "supplemental") {
        """
        Rscript ${projectDir}/bin/wgs/annotate-bedpe-with-genes.r \
            --ensembl=${params.ensemblUniqueBed} \
            --bedpe=${annot_sv_bedpe} \
            --out_file=${sampleID}.MDLS_sv_annotated_genes_supplemental.bed \
            --supplemental
        """
    }
}
