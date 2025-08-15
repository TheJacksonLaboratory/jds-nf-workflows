process SV_MERGE {
    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.bedpe", mode:'copy'

    input:
        tuple val(sampleID), file(vcf_tuple)
        val(chrom_list)

    output:
        tuple val(sampleID), file("*.mergedCall.DLMS.bedpe"), emit: bedpe
        tuple val(sampleID), file("*.mergedCall.DLMS.supplemental.bedpe"), emit: supp_bedpe

    script:
    
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    """
        Rscript ${projectDir}/bin/germline_sv/merge_sv.r \
        --vcf=${vcf_tuple[0]},${vcf_tuple[1]},${vcf_tuple[2]},${vcf_tuple[3]} \
        --caller=delly,lumpy,manta,svaba \
        --sample_name=${sampleID} \
        --build=${params.genome_build} \
        --slop=1000 \
        --allowed_chr=${listOfChroms} \
        --min_sv_length=200 \
        --out_file=${sampleID}.mergedCall.DLMS.bedpe \
        --out_file_supplemental=${sampleID}.mergedCall.DLMS.supplemental.bedpe
    """
}
