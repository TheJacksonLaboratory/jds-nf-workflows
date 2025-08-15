process BCFTOOLS_QUERY_DELLY_CNV {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.bed", mode: 'copy'

    input:
    tuple val(sampleID), path(bcf), path(csi)

    output:
    tuple val(sampleID), path("*segmentation.bed"), emit: segmentation_file

    script:

    """
    bcftools query -H -f "%CHROM\\t%POS\\t%INFO/END\\t%ID\\t[%RDCN]\\t[%RDSD]\\t[%CN]\n" ${bcf} > ${sampleID}_delly_cnv_segmentation.bed
    """
}

// WARN: Input tuple does not match tuple declaration in process `WGS:WGS_SV:BCFTOOLS_QUERY_DELLY_CNV` -- offending value: [hg38_WGS, /projects/compsci/omics_share/meta/benchmarking/wgs_sv/work/7d/70a1b2f02fa89a2b6773b9a201208e/hg38_WGS_delly_cnv_classified.bcf, /projects/compsci/omics_share/meta/benchmarking/wgs_sv/work/7d/70a1b2f02fa89a2b6773b9a201208e/hg38_WGS_delly_cnv_classified.bcf.csi]