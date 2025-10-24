process SV_MERGE {
    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID}/merged_sv", pattern:"*.bedpe", mode:'copy'

    input:
        tuple val(sampleID), file(vcf_tuple)
        val(chrom_list)

    output:
        tuple val(sampleID), file("*.mergedCall.MDLS.bedpe"), emit: bedpe
        tuple val(sampleID), file("*.mergedCall.MDLS.supplemental.bedpe"), emit: supp_bedpe

    script:
    
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    if (params.workflow == "wgs" || params.workflow == "wgs_sv_bam") {

        """
            Rscript ${projectDir}/bin/wgs/merge_sv.r \
            --vcf=${vcf_tuple[0]},${vcf_tuple[1]},${vcf_tuple[2]},${vcf_tuple[3]} \
            --callers=manta,delly,lumpy,svaba \
            --sample_name=${sampleID} \
            --build=${params.genome_build} \
            --slop=${params.sv_slop} \
            --sizemargin=${params.sizemargin} \
            --allowed_chr=${listOfChroms} \
            --min_sv_length=${params.min_sv_length} \
            --out_file=${sampleID}.mergedCall.MDLS.bedpe \
            --out_file_supplemental=${sampleID}.mergedCall.MDLS.supplemental.bedpe
        """

    } else if (params.workflow == "wgs_long_read") {

        """
            Rscript ${projectDir}/bin/wgs/merge_sv.r \
            --vcf=${vcf_tuple[0]},${vcf_tuple[1]} \
            --callers=sniffles,pbsv \
            --sample_name=${sampleID} \
            --build=${params.genome_build} \
            --slop=${params.sv_slop} \
            --sizemargin=${params.sizemargin} \
            --allowed_chr=${listOfChroms} \
            --min_sv_length=${params.min_sv_length} \
            --out_file=${sampleID}.mergedCall.MDLS.bedpe \
            --out_file_supplemental=${sampleID}.mergedCall.MDLS.supplemental.bedpe
        """
    
    }
}
