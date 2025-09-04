process ANNOTATE_SV_WITH_CNV {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.bedpe", mode: 'copy'
    
    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    input:
        tuple val(sampleID), val(normal_name), val(tumor_name), file(delly_annot), file(annot_sv_genes_bedpe)
        val(suppl_switch)

    output:
        tuple val(sampleID), file("*.bedpe"), val(normal_name), val(tumor_name), emit: sv_genes_cnv_bedpe
    
    script:
        if (suppl_switch == "main") {
            output_name = "${sampleID}_manta_lumpy_delly_svaba_sv_annotated_genes_cnv" == annot_sv_genes_bedpe.baseName ? "${sampleID}_manta_lumpy_delly_svaba_sv_reannotated_genes_cnv.bedpe" : "${sampleID}_manta_lumpy_delly_svaba_sv_annotated_genes_cnv.bedpe"
            """
            Rscript ${projectDir}/bin/pta/annotate-bedpe-with-cnv.r \
                --cnv=${delly_annot} \
                --bedpe=${annot_sv_genes_bedpe} \
                --distance_limit=1000 \
                --out_file=${output_name}
            """
        } else if (suppl_switch == "supplemental") {
            output_name = "${sampleID}_manta_lumpy_delly_svaba_sv_annotated_genes_cnv_supplemental" == annot_sv_genes_bedpe.baseName ? "${sampleID}_manta_lumpy_delly_svaba_sv_reannotated_genes_cnv_supplemental.bedpe" : "${sampleID}_manta_lumpy_delly_svaba_sv_annotated_genes_cnv_supplemental.bedpe"
            """
            Rscript ${projectDir}/bin/pta/annotate-bedpe-with-cnv.r \
                --cnv=${delly_annot} \
                --bedpe=${annot_sv_genes_bedpe} \
                --distance_limit=1000 \
                --out_file=${output_name}
            """
        }
}
