process MITY_RUNALL {
    tag "${sampleID}"

    cpus 1
    memory 40.GB
    time "2:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/mity:2.0.0--pyhdfd78af_0'

    publishDir "${params.pubdir}/${sampleID}/mt_callers", pattern: "*.vcf.gz*", mode:'copy'
    publishDir "${params.pubdir}/${sampleID}/mt_callers", pattern: "*.xlsx", mode:'copy', enabled: params.gen_org == 'human'

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("*.normalise.vcf.gz"), path("*.normalise.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(sampleID), path("*.mity.report.xlsx"), emit: report

    script:
    reference = params.gen_org == 'mouse' ? 'mm10' : 'hg38'
    """
    mity runall \
    --prefix ${sampleID}.mity \
    --reference ${reference} \
    --min_vaf 0.01 \
    --contig ${params.mt_contig_name} \
    --custom-reference-fasta ${params.mt_fasta} \
    --custom-reference-genome ${params.mt_genome} \
    ${bam}
    """
}
