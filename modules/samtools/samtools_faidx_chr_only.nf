process SAMTOOLS_FAIDX_CHR_ONLY {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    //publishDir "${params.pubdir}/g2gtools", pattern:"*.fai", mode:'copy', enabled: params.workflow == 'generate_pseudoreference' ? true : false
    //publishDir "${params.pubdir}/genome_info", pattern:"*.fai", mode:'copy', enabled: params.workflow == 'chipseq' ? true : false

    input:
        path(fasta)

    output:
        path("*CHR_Only.fa"), emit: chr_fa

    script:
    suffix = params.gen_org=='mouse' ? "\$(echo -e ${fasta} | sed 's/.toplevel.fa//g')" : "\$(echo -e ${fasta} | sed 's/_assembly38.fasta//g')"

    if (params.gen_org=='mouse')
    """
      samtools faidx ${fasta} 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT > ${suffix}.primary_CHR_Only.fa
    """
    else if (params.gen_org=='human')
    """
      samtools faidx ${fasta} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrM chrX chrY > ${suffix}.${params.genome_build}.dna.CHR_Only.fa
    """
}
