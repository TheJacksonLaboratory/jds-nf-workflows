process ANNOTATE_DELLY_CNV {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.bed", mode: 'copy'

    input:
        tuple val(sampleID), file(delly_cnv)
        val(chrom_list)

    output:
        tuple val(sampleID), file("${sampleID}_cnv_annotated_final.bed"), emit: delly_annot
        tuple val(sampleID), file("${sampleID}_cnv_annotated_supplemental.bed"), emit: delly_annot_suppl

    script:
        listOfChroms = chrom_list.collect { "$it" }.join(',')

        if ( params.gen_org == 'mouse' )
        """
        Rscript ${projectDir}/bin/wgs/annotate-cnv-delly.r \
            --cnv=${delly_cnv} \
            --caller="delly" \
            --sample_name=${sampleID} \
            --cytoband=${params.cytoband} \
            --db_names="INS,DEL,INV" \
            --db_files=${params.known_del},${params.known_ins},${params.known_inv} \
            --ensembl=${params.ensemblUniqueBed} \
            --allowed_chr=${listOfChroms} \
            --overlap_fraction=0.8 \
            --out_file_main=${sampleID}_cnv_annotated_final.bed \
            --out_file_supplemental=${sampleID}_cnv_annotated_supplemental.bed
        """

        else if ( params.gen_org == 'human' )
        """
        Rscript ${projectDir}/bin/wgs/annotate-cnv-delly.r \
            --cnv=${delly_cnv} \
            --caller="delly" \
            --sample_name=${sampleID} \
            --cytoband=${params.cytoband} \
            --db_names="DGV,1000G,COSMIC" \
            --db_files=${params.dgv},${params.thousandG},${params.cosmicUniqueBed} \
            --ensembl=${params.ensemblUniqueBed} \
            --allowed_chr=${listOfChroms} \
            --overlap_fraction=0.8 \
            --out_file_main=${sampleID}_cnv_annotated_final.bed \
            --out_file_supplemental=${sampleID}_cnv_annotated_supplemental.bed

        """
}


